#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

__all__ = []

import kplr
import emcee
import transit
import numpy as np
from functools import partial
import matplotlib.pyplot as pl
from emcee.autocorr import integrated_time

import george
from george import kernels


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
        m = np.array(m * (qual[chunk] == 0), dtype=bool)

        # Skip it if there aren't enough data points.
        if sum(m) <= min_length:
            continue

        # Normalize the flux to have a median of 1.
        mu = np.median(flux[chunk][m])
        light_curves.append(map(np.ascontiguousarray,
                                (time[chunk][m], flux[chunk][m] / mu,
                                 ferr[chunk][m] / mu)))

    return light_curves


class ProbabilisticModel(object):

    def __init__(self, system, planet, gp_models):
        self.system = system
        self.planet = planet
        self.gp_models = gp_models

    def get_parameters(self):
        return np.array([
            np.log(self.system.central.radius),
            self.system.central.q1,
            self.system.central.q2,
            np.log(self.planet.period),
            self.planet.t0,
            np.log(self.planet.r),
            self.planet.b,
        ])

    def set_parameters(self, p):
        # Update the stellar parameters.
        self.system.central.radius = np.exp(p[0])
        self.system.central.q1 = p[1]
        self.system.central.q2 = p[2]

        # And the planet parameters.
        self.planet.period = np.exp(p[3])
        self.planet.t0 = p[4]
        self.planet.r = np.exp(p[5])
        self.planet.b = p[6]

    def lnprior(self, p):
        if not (0 < p[1] < 1 and 0 < p[2] < 1):  # LD
            return -np.inf

        if not 0 < p[6] < 1.5:  # b
            return -np.inf

        # Radius measurement.
        return -0.5 * (np.exp(p[0]) - 0.95) ** 2 / 0.1 ** 2

    def lnlike(self):
        return sum(gp.lnlikelihood(y, quiet=True) for gp, y in models)

    def lnprob(self, p):
        lp = self.lnprior(p)
        if not np.isfinite(lp):
            return -np.inf
        try:
            self.set_parameters(p)
        except ValueError:
            return -np.inf
        return self.lnlike() + lp

    def __call__(self, p):
        return self.lnprob(p)


if __name__ == "__main__":
    # Download the data.
    client = kplr.API()
    kicid = 8644545
    kic = client.star(kicid)
    datasets = kic.get_light_curves()

    # Process the light curves.
    light_curves = [lc for lcs in map(prepare_light_curve, datasets)
                    for lc in lcs]

    # Throw away the datasets without transits.
    period, t0 = 295.963, 138.91
    hp = 0.5 * period
    selection = lambda lc: np.any(np.abs((lc[0]-t0+hp) % period-hp) < 1.0)
    light_curves = filter(selection, light_curves)
    print(len(light_curves))

    # Plot the raw data.
    for lc in light_curves:
        pl.plot(lc[0], lc[1], ".", ms=3)
    pl.savefig("raw_data.png")

    # Set up the initial system.
    system = transit.System(transit.Central(radius=0.95))
    planet = transit.Body(r=2.03 * 0.01, period=period, t0=t0, b=0.9)
    system.add_body(planet)
    texp = kplr.EXPOSURE_TIMES[1] / 60. / 60. / 24.
    mean_function = partial(system.light_curve, texp=texp)

    # Set up the Gaussian processes.
    pl.clf()
    offset = 0.001
    models = []
    for i, lc in enumerate(light_curves):
        dt = np.median(np.diff(lc[0])) * integrated_time(lc[1])
        kernel = np.var(lc[1]) * kernels.Matern32Kernel(dt ** 2)
        gp = george.GP(kernel, mean=mean_function, solver=george.HODLRSolver)
        gp.compute(lc[0], lc[2])
        models.append((gp, lc[1]))

        t = (lc[0]-t0+hp) % period-hp
        pl.plot(t, lc[1] + i * offset, ".k", ms=3)
        pl.plot(t, gp.predict(lc[1], lc[0], mean_only=True) + i*offset, "b")

    pl.savefig("initial.png")
    pl.xlim(-5, 5)
    pl.savefig("initial_zoom.png")

    model = ProbabilisticModel(system, planet, models)
    p0 = model.get_parameters()

    ndim = len(p0)
    nwalkers = 32
    p0 = [p0 + 1e-6 * np.random.randn(ndim) for i in range(nwalkers)]

    sampler = emcee.EnsembleSampler(nwalkers, ndim, model)
    print("Running burn-in")
    p0, _, _ = sampler.run_mcmc(p0, 1000)
    sampler.reset()
    print("Running production")
    sampler.run_mcmc(p0, 1000)

    for i in range(ndim):
        pl.clf()
        pl.plot(sampler.chain[:, :, i].T, color="k", alpha=0.3)
        pl.savefig("samples/{0:03d}.png".format(i+1))
