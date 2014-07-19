.. _quickstart:

Getting started
===============

Installation
------------

You can install the most recent stable version of transit using `PyPI
<#stable>`_ or the development version from `GitHub <#dev>`_.

Prerequisites
+++++++++++++

Whichever method you choose, you'll need to make sure that you first have
`Boost <http://www.boost.org/>`_ installed.
On Linux:

.. code-block:: bash

    sudo apt-get install libboost-dev

On Mac:

.. code-block:: bash

    brew install boost

.. note:: Chances are high that transit won't work on Windows right now
    because  it hasn't been tested at all but feel free to try it out at your
    own risk!

You'll also need a working scientific Python installation (including `NumPy
<http://www.numpy.org/>`_).
I recommend the `Anaconda distribution <http://continuum.io/downloads>`_ if
you don't already have your own opinions.

.. _stable:

Stable Version
++++++++++++++

The simplest way to install the `most recent stable version of transit
<https://pypi.python.org/pypi/trasit>`_ is using `pip
<http://pip.readthedocs.org/>`_:

.. code-block:: bash

    pip install transit

If you installed Boost in a strange place, specify that location by running
(sorry to say that it's pretty freaking ugly):

.. code-block:: bash

    pip install transit \
        --global-option=build_ext \
        --global-option=-I/path/to/boost


.. _dev:

Development Version
+++++++++++++++++++

To get the source for the development version, clone the git repository:

.. code-block:: bash

    git clone https://github.com/dfm/transit.git
    cd transit

Then, install the package by running the following command:

.. code-block:: bash

    python setup.py install

If installed Boost in a non-standard location, you can specify the correct
path using the install command:

.. code-block:: bash

    python setup.py build_ext -I/path/to/boost install

Testing
+++++++

To run the unit tests, install `nose <https://nose.readthedocs.org>`_ and then
execute:

.. code-block:: bash

    nosetests transt -v

All of the tests should (of course) pass.
If any of the tests don't pass and if you can't sort out why, `open an issue
on GitHub <https://github.com/dfm/transit/issues>`_.


A Simple Example
----------------

In this example, we'll build a simple Earth-like transit and plot the result:

.. code-block:: python

    import transit
    import numpy as np
    import matplotlib.pyplot as pl

    # Build the transiting system.
    s = transit.System(transit.Central())
    body = transit.Body(r=0.009155, period=365.25, t0=0.99, b=0.2, e=0.0167)
    s.add_body(body)

    # Compute the light curve integrated over a Kepler long cadence
    # exposure time.
    texp = 1626.0 / 86400.0
    t = np.arange(0, 2, texp)
    f = s.light_curve(t, texp=texp)

    # Plot the results
    pl.plot(t, (f-1) * 1e6, ".k")

This should produce a figure that looks something like this:

.. image:: _static/simple.png

For a list of all the available options, see :ref:`api`.
