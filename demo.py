#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import time
import numpy as np
import matplotlib.pyplot as pl

from kplr import EXPOSURE_TIMES

import transit

texp = EXPOSURE_TIMES[1] / 86400.0

s = transit.System(transit.Central())
body = transit.Body(r=0.02, mass=9e-4, period=3.0, t0=5, b=0.0, e=0.3,
                    pomega=0.5)
s.add_body(body)

t = np.linspace(0, 10.0, 1000)

fig, axes = pl.subplots(2, 1)
solver = s._get_solver()

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
# print(time.time() - strt)
axes[0].plot(t, f, ".k")
axes[0].axvline(body.t0, color="k")

# strt = time.time()
# f = s.light_curve(t, texp=0.2)
# print(time.time() - strt)
# axes[0].plot(t, f, ".r")

axes[1].plot(t, s.radial_velocity(t), ".k")
axes[1].axvline(body.t0, color="k")

fig.savefig("demo.pdf")
