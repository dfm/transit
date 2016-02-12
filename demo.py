#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import time
import numpy as np
import matplotlib.pyplot as pl

import emcee

from kplr import EXPOSURE_TIMES

import transit

texp = EXPOSURE_TIMES[1] / 86400.0

s = transit.System(transit.Central(dilution=0.05))
body = transit.Body(radius=0.2, mass=0.0, period=4.0, t0=2, b=0.0, e=0.4,
                    omega=0.5*np.pi + 0.01)
s.add_body(body)

s.thaw_parameter("*")
print(s.get_parameter_names())

x = np.linspace(0, 10.0, 1000)
yerr = 5e-4 * np.ones_like(x)
y = s.light_curve(x) + yerr * np.random.randn(len(x))

pl.plot(x, y, ".k")
pl.savefig("data.png")

p0 = s.get_vector()

assert 0


# fig, axes = pl.subplots(2, 1)

# strt = time.time()
# f1 = s.light_curve(t)
# print(time.time() - strt)
# pl.plot(t, f1, "k")

# s.central.dilution = 0.5
# body.radius *= 1.4
# f2 = s.light_curve(t)
# # f2 = (f2 - 1.0) * 2 + 1
# print(f1 - f2)
# pl.plot(t, f2, "b")

# pl.gca().axvline(body.t0, color="k")
# pl.savefig("face.png")
# assert 0

# solver = s._get_solver()
# assert 0

eps = 1e-7
p = solver.position(t+eps)
m = solver.position(t-eps)
vel = solver.velocity(t)
# print(vel)
# print(0.5 * (p - m) / eps)
print(np.abs(vel - 0.5 * (p - m) / eps).max())
# assert 0

# print((0.5 * (p - m) / eps, vel))

# assert 0

strt = time.time()
f = s.light_curve(t)
print(time.time() - strt)
axes[0].plot(t, f, "k")
axes[0].axvline(body.t0, color="k")

# strt = time.time()
# f = s.light_curve(t, texp=0.2)
# print(time.time() - strt)
# axes[0].plot(t, f, ".r")

axes[1].plot(t, s.radial_velocity(t), "k")
axes[1].axvline(body.t0, color="k")
axes[1].axhline(0, color="k")

fig.savefig("demo.pdf")
