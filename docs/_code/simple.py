#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as pl

import transit

# Build the transiting system.
s = transit.System(transit.Central())
body = transit.Body(r=0.009155, period=365.25, t0=0.99, b=0.2, e=0.0167)
s.add_body(body)

# Compute the light curve integrated over a Kepler long cadence exposure time.
texp = 1626.0 / 86400.0
t = np.arange(0, 2, texp)
f = s.light_curve(t, texp=texp)

# Plot the result.
fig = pl.figure(figsize=(6, 3.5))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.2, bottom=0.2, right=0.95, top=0.95)
ax.plot(t, (f-1) * 1e6, ".k")

ax.set_xlim(0, 2)
ax.set_ylim(-118, 8)
ax.set_ylabel("relative flux [ppm]")
ax.set_xlabel("time [days]")

fig.savefig("../_static/simple.png", dpi=150)
