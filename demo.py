#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as pl

import transit

s = transit.System(transit.Central())
s.add_body(transit.Body(flux=0.01, r=0.02, period=10.0, t0=3.0))

t = np.linspace(0, 50, 5000)
f = s.light_curve(t)

pl.plot(t, f)
pl.savefig("demo.png")
