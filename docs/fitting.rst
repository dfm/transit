.. module:: transit

.. _api:

Fitting a Kepler Light Curve
============================

This is an example where we'll be fitting for the parameters of a transiting
exoplanet using real Kepler data.

.. note:: **Prerequisites** — You'll need `fitsio
    <https://github.com/esheldon/fitsio>`_, `kplr <http://dfm.io/kplr>`_,
    `George <http://dfm.io/george>`_, and `emcee <http://dfm.io/emcee>`_
    installed.

Getting the data
----------------

Let's fit a small planet on a long orbit from `Petigura's sample
<http://arXiv.org/abs/1311.6806>`_:

.. code-block:: python

    import kplr
    client = kplr.API()

    kicid = 8644545
    kic = client.star(kicid)
    datasets = kic.get_light_curves()

Then, we'll write a function that takes each of these light curves and break
them up into chunks that each approximately one month long (breaking at gaps
in the data):

.. code-block:: python

    import numpy as np

    def prepare_light_curve(lc, tol=20, min_length=0):
        data = lc.read(columns=["TIME", "SAP_FLUX", "SAP_FLUX_ERR", "SAP_QUALITY"])
        time = data["TIME"]
        flux = data["SAP_FLUX"]
        ferr = data["SAP_FLUX_ERR"]
        qual = data["SAP_QUALITY"]

        # Loop over the time array and break it into "chunks" when there is "a
        # sufficiently long gap" with no data.
        count, current, chunks = 0, [], []
        for i, t in enumerate(time):
            if np.isnan(t):
                count += 1
            else:
                if count > tol:
                    chunks.append(list(current))
                    current = []
                    count = 0
                current.append(i)
        if len(current):
            chunks.append(current)

        # Loop over the chunks and construct the output.
        light_curves = []
        for chunk in chunks:
            # Remove missing or bad data points.
            m = np.prod(map(np.isfinite, (time[chunk], flux[chunk], ferr[chunk])),
                        axis=0)
            m = np.array(m * (qual[chunk == 0]), dtype=bool)

            # Skip it if there aren't enough data points.
            if sum(m) <= min_length:
                continue

            # Normalize the flux to have a median of 1.
            mu = np.median(flux[chunk][m])
            light_curves.append(map(np.ascontiguousarray,
                                    (time[chunk][m], flux[chunk][m] / mu,
                                     ferr[chunk][m] / mu)))

        return light_curves

And apply the function to the light curves that we downloaded:

.. code-block:: python

    light_curves = [lc for lcs in map(prepare_light_curve, datasets) for lc in lcs]

This planet has a long period (296 days) so we'll conservatively throw away
datasets that probably don't have transits:

.. code-block:: python

    period, t0 = 295.963, 138.91
    hp = 0.5 * period
    selection = lambda lc: np.any(np.abs((lc[0]-t0+hp) % period-hp) < 2.0)
    light_curves = filter(selection, light_curves)

Noise model
-----------

We'll use a Gaussian process noise model for each of these light curves and
marginalize over the hyperparameters.

..
