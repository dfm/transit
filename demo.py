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
body = transit.Body(flux=5e-5, r=0.02, period=100.0, t0=0.99, b=0.2, e=0.2)
s.add_body(body)

t = np.linspace(0, 2, 1000)

strt = time.time()
f = s.light_curve(t)
print(time.time() - strt)
pl.plot(t, f)

strt = time.time()
f = s.light_curve(t, texp=0.2)
print(time.time() - strt)
pl.plot(t, f)

pl.savefig("demo.pdf")
