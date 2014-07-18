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
s.add_body(transit.Body(flux=5e-5, r=0.02, period=2.0, t0=0.0))

t = np.linspace(0, 4, 1000)
strt = time.time()
f = s.light_curve(t, texp=texp)
print(time.time() - strt)

pl.plot(t, f)
pl.savefig("demo.pdf")
