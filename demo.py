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
body = transit.Body(r=0.02, mass=9e-4, period=1.0, t0=0.5, b=0.2, e=0.9,
                    pomega=0.0)
s.add_body(body)

t = np.linspace(0, 1.0, 1000)

fig, axes = pl.subplots(2, 1)
solver = s._get_solver()

eps = 1e-8
p = solver.position(t+eps)
m = solver.position(t-eps)
vel = solver.velocity(t)
print(vel - 0.5 * (p - m) / eps)
assert 0

# print((0.5 * (p - m) / eps, vel))

# assert 0

strt = time.time()
f = s.light_curve(t)
print(time.time() - strt)
axes[0].plot(t, f, ".k")

# strt = time.time()
# f = s.light_curve(t, texp=0.2)
# print(time.time() - strt)
# axes[0].plot(t, f, ".r")

# axes[1].plot(t, s.radial_velocity(t), ".k")

fig.savefig("demo.pdf")
