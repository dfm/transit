# -*- coding: utf-8 -*-

from __future__ import division, print_function

__all__ = ["SimpleSystem"]

import numpy as np
from ._transit import CythonSolver


class SimpleSystem(object):

    def __init__(self, period=None, ror=None, duration=None, t0=0.0,
                 impact=0.0, flux=1.0, q1=None, q2=None, mu1=None, mu2=None):
        if period is None or ror is None or duration is None:
            raise ValueError("missing required parameter 'period', 'ror', or "
                             "'duration'")
        self.period = period
        self.t0 = t0
        self.duration = duration
        self.ror = ror
        self.impact = impact
        self.flux = flux
        self.unfrozen = np.ones(7, dtype=bool)

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
        return CythonSolver().simple_light_curve(self._get_params(),
                                                 t, texp, tol, maxdepth)

    def light_curve_gradient(self, t, texp=0.0, tol=1e-8, maxdepth=4):
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
        f, df = CythonSolver().simple_gradient(self._get_params(),
                                               t, texp, tol, maxdepth)

        return f, df

    def __len__(self):
        return np.sum(self.unfrozen)

    def _parameter_names(self):
        return ["ln_ror", "ln_period", "t0", "impact", "ln_duration",
                "q1_param", "q2_param"]

    def get_parameter_names(self, full=False):
        if full:
            return self._parameter_names()
        return [n for n, f in zip(self._parameter_names(), self.unfrozen)
                if f]

    def _get_params(self):
        return np.array([
            np.log(self.ror), np.log(self.period), self.t0, self.impact,
            np.log(self.duration),
            np.log(self.q1)-np.log(1.0-self.q1),
            np.log(self.q2)-np.log(1.0-self.q2),
        ])

    def get_vector(self):
        return self._get_params()[self.unfrozen]

    def check_vector(self, *args):
        return True

    def set_vector(self, vector):
        p = self._get_params()
        p[self.unfrozen] = vector
        self.ror = np.exp(p[0])
        self.period = np.exp(p[1])
        self.t0 = p[2]
        self.impact = p[3]
        self.duration = np.exp(p[4])
        self.q1 = max(0.0, min(1.0, 1.0 / (1. + np.exp(-p[5]))))
        self.q2 = max(0.0, min(1.0, 1.0 / (1. + np.exp(-p[6]))))

    def get_value(self, t, **kwargs):
        return self.light_curve(t, **kwargs)

    def get_gradient(self, t, **kwargs):
        return self.light_curve_gradient(t, **kwargs)[1][:, self.unfrozen].T

    def get_parameter(self, *args):
        raise NotImplementedError()

    def set_parameter(self, *args):
        raise NotImplementedError()

    def freeze_parameter(self, parameter_name):
        i = self._parameter_names().index(parameter_name)
        self.unfrozen[i] = False

    def thaw_parameter(self, parameter_name):
        i = self._parameter_names().index(parameter_name)
        self.unfrozen[i] = True

    def get_bounds(self):
        return [(None, None) for _ in range(len(self))]
