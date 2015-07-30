# -*- coding: utf-8 -*-

from __future__ import division, print_function

__all__ = ["SimpleSystem"]

import numpy as np
from ._transit import PythonSimpleSolver as Solver


class SimpleSystem(object):

    def __init__(self, period=None, ror=None, duration=None, t0=0.0,
                 impact=0.0, flux=1.0,
                 q1=None, q2=None, mu1=None, mu2=None):
        if period is None or ror is None or duration is None:
            raise ValueError("missing required parameter 'period', 'ror', or "
                             "'duration'")
        self.period = period
        self.t0 = t0
        self.duration = duration
        self.ror = ror
        self.impact = impact
        self.flux = flux

        # Allow different limb darkening parameters.
        if mu1 is not None and mu2 is not None:
            if q1 is not None or q2 is not None:
                raise RuntimeError("You can't use *both* limb-darkening "
                                   "parameterizations!")
            self.coeffs = (mu1, mu2)
        else:
            self.q1 = q1 if q1 is not None else 0.5
            self.q2 = q2 if q2 is not None else 0.5

    @property
    def q1(self):
        return self._q1

    @q1.setter
    def q1(self, v):
        if not 0 <= v <= 1:
            raise ValueError("Invalid limb darkening coefficient")
        self._q1 = v

    @property
    def q2(self):
        return self._q2

    @q2.setter
    def q2(self, v):
        if not 0 <= v <= 1:
            raise ValueError("Invalid limb darkening coefficient")
        self._q2 = v

    @property
    def coeffs(self):
        q1, q2 = self.q1, self.q2
        q1 = np.sqrt(np.abs(q1))
        return 2*q1*q2, q1*(1-2*q2)

    @coeffs.setter
    def coeffs(self, value):
        u1, u2 = value
        u2 = u1+u2
        self.q1, self.q2 = u2*u2, 0.5*u1/u2

    def _get_solver(self):
        u1, u2 = self.coeffs
        return Solver(u1, u2, self.period, self.t0, self.duration, self.ror,
                      self.impact)

    def light_curve(self, t, texp=0.0, tol=1e-8, maxdepth=4):
        """
        Get the light curve evaluated at a list of times using the current
        model.

        :param t:
            The times where the light curve should be evaluated (in days).

        :param tol:
            The stopping criterion for the exposure time integration.

        :param maxdepth:
            The maximum recursion depth of the exposure time integrator.

        """
        t = np.atleast_1d(t)
        return self._get_solver().light_curve(self.flux, t, texp, tol,
                                              maxdepth)
