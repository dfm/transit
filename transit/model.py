# -*- coding: utf-8 -*-

from __future__ import division, print_function

__all__ = ["Central", "Body", "System"]

try:
    from itertools import izip, imap
except ImportError:
    izip, imap = zip, map

import numpy as np

from ._transit import ldlc_kepler


# Newton's constant in $R_\odot^3 M_\odot^{-1} {days}^{-2}$.
_G = 2945.4625385377644


class Central(object):
    """
    The "central"---in this context---is the massive central body in a
    :class:`System`.

    :param mass:
        The mass of the body measured in Solar masses. (default: ``1.0``)

    :param radius:
        The radius of the body measured in Solar radii. (default: ``1.0``)

    :param flux:
        The un-occulted flux measured in whatever units you feel like using.
        (default: ``1.0``)

    **Limb darkening** can be specified using ``(mu1, mu2)`` or ``(q1, q2)``.
    TODO: explain.

    """

    def __init__(self, mass=1.0, radius=1.0, flux=1.0, q1=None, q2=None,
                 mu1=None, mu2=None):
        self.mass = mass
        self.radius = radius
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
    def coeffs(self):
        q1, q2 = self.q1, self.q2
        q1 = np.sqrt(np.abs(q1))
        return 2*q1*q2, q1*(1-2*q2)

    @coeffs.setter
    def coeffs(self, value):
        u1, u2 = value
        u2 = u1+u2
        self.q1, self.q2 = u2*u2, 0.5*u1/u2

    @property
    def density(self):
        r = self.radius
        return 0.75 * self.mass / (np.pi * r * r * r)

    @density.setter
    def density(self, rho):
        r = self.radius
        self.mass = np.pi * rho * r * r * r / 0.75


class Body(object):
    r"""
    A "body"---in this context---is a (possibly) massive body orbiting a
    :class:`Central` in a :class:`System`. There are several ways to
    initialize this and once it has been added to a system using the
    :func:`System.add_body` method, they should all be equivalent. The orbital
    elements specified either specify a Keplerian orbit. This object includes
    all sorts of magic for converting between different specifications when
    needed but the base description of the planet and the orbit is
    parameterized by the parameters:

    .. code-block:: python

        (flux, r, mass, a, t0, e, pomega, ix, iy)

    :param flux:
        The flux of the body measured relative to the central.
        (default: ``0.0``)

    :param r:
        The radius measured in Solar radii. (default: ``0.0``)

    :param mass:
        The mass in Solar masses. (default: ``0.0``)

    :param a:
        The semi-major axis of the orbit measured in Solar radii. Either this
        parameter or ``period`` must be provided but not both.

    :param period:
        The period of the orbit in days. Either this parameter or ``a`` must
        be provided but not both.

    :param t0:
        The epoch of the orbit in days. (default: ``0.0``)

    :param e:
        The eccentricity of the orbit. (default: ``0.0``)

    :param pomega:
        The orientation of the orbital ellipse in radians. (default: ``0.0``)

    :param ix:
        The relative inclination of the orbital plane along the line-of-sight
        in degrees. This angle is measured differently than you're used to:
        zero degrees is edge on and 90 degrees in face on. This angle will be
        subtracted from the base inclination of the planetary system to get
        the standard measurement of the inclination. Either this parameter
        or ``b`` can be specified but not both. (default: ``0.0``)

    :param b:
        The mean impact parameter of the orbit measured in stellar radii (not
        Solar radii). Specifically, this impact parameter is defined as

        .. math::

            b = \frac{a}{R_\star} \tan i_x

        (default: ``0.0``)

    :param iy:
        The rotation of the orbital ellipse in the plane of the sky measured
        in radians. Note: this value will not affect the light curve or radial
        velocity values at all. (default: ``0.0``)

    """

    def __init__(self, flux=0.0, r=0.0, mass=0.0, a=None, period=None, t0=0.0,
                 e=0.0, pomega=0.0, ix=None, b=None, iy=0.0):
        self.flux = flux
        self.r = r
        self._a = a
        self._period = period
        self.mass = mass
        self.t0 = t0
        self.e = e
        self.pomega = pomega
        self._b = b
        self._ix = ix
        self.iy = iy

        if a is not None:
            assert period is None, \
                "You can't supply both a period and a semi-major axis."
        else:
            assert period is not None, \
                "You must supply either a period or a semi-major axis."

        assert b is None or ix is None, \
            "You can't supply both an inclination and impact parameter."
        if ix is None and b is None:
            self._ix = 0.0

    def _check_ps(self):
        if not hasattr(self, "system"):
            raise RuntimeError("You must add this body to a system "
                               "before getting the period.")

    @property
    def period(self):
        # If we already have a period, return that.
        if self._period is not None:
            return self._period

        # If not, check to make sure that we're already part of a system
        # and then compute the period based on the star's mass.
        self._check_ps()
        mstar = self.system.central.mass
        a = self._a
        return 2 * np.pi * np.sqrt(a * a * a / _G / (mstar + self.mass))

    @period.setter
    def period(self, P):
        self._check_ps()
        mstar = self.system.central.mass
        self._a = (_G*P*P*(self.mass+mstar)/(4*np.pi*np.pi)) ** (1./3)
        self._period = None

    @property
    def a(self):
        if self._a is None:
            self.period = self._period
        return self._a

    @a.setter
    def a(self, a):
        self._period = None
        self._a = a

    @property
    def b(self):
        # If we already have an impact parameter, return that.
        if self._b is not None:
            return self._b

        # If not, check to make sure that we're already part of a system
        # and then compute the impact parameter based on the star's radius.
        self._check_ps()
        rstar = self.system.central.radius
        return self.a * np.tan(np.radians(self._ix)) / rstar

    @b.setter
    def b(self, b):
        self._check_ps()
        rstar = self.system.central.radius
        self._ix = np.degrees(np.arctan2(b, self.a / rstar))
        self._b = None

    @property
    def ix(self):
        if self._ix is None:
            self.b = self._b
        return self._ix

    @ix.setter
    def ix(self, ix):
        self._b = None
        self._ix = ix


class System(object):
    """
    A "planetary system" contains a "central" bright :class:`Central` and some
    number (``>= 0``) :class:`Body` objects orbiting. The default orbits
    are purely Keplerian but sub-classes can include more sophisticated
    solvers.

    :param central:
        A :class:`Central` that specifies the central bright object.

    :para iobs:
        The inclination of the mean orbital plane in degrees. This is
        measured in the standard way with zero inclination meaning face on and
        90 degrees is edge on. (default: ``90.0``)

    """

    def __init__(self, central, iobs=90.0):
        self.central = central
        self.central.system = self
        self.bodies = []
        self.iobs = iobs

    def __len__(self):
        return len(self.bodies)

    def add_body(self, body):
        """
        Add a :class:`Body` to the system. This function also sets the
        ``system`` attribute of the body.

        :param body:
            The :class:`Body` to add.

        """
        body.system = self
        self.bodies.append(body)

    def _get_pars(self):
        """
        Coerce the planetary parameters into the form needed by the
        :func:`lightcurve` function. This function returns a list of the
        following arrays

        .. code-block:: python

            flux, mass, r, a, t0, e, pomega, ix, iy

        """
        r = ((p.flux, p.mass, p.r, p.a, p.t0, p.e, p.pomega, p.ix, p.iy)
             for p in self.bodies)
        return izip(*r)

    def light_curve(self, t, texp=0.0, tol=0.0, maxdepth=0):
        """
        Get the light curve evaluated at a list of times using the current
        model.

        :param t:
            The times where the light curve should be evaluated (in days).

        """
        if len(self) == 0:
            return self.central.flux + np.zeros_like(t)

        # Coerce the planetary orbital parameters.
        pars = list(imap(np.atleast_1d, self._get_pars()))
        flux, mass, r, a, t0, e, pomega, ix, iy = pars

        # Grab some shortcuts to some useful objects.
        s = self.central
        t = np.atleast_1d(t)
        mu1, mu2 = s.coeffs

        # Compute the light curve using the Kepler solver.
        lc = ldlc_kepler(t, mu1, mu2, s.mass, s.radius, flux, mass, r, a, t0,
                         e, pomega, np.radians(90.-self.iobs+ix),
                         np.radians(iy), texp, tol, maxdepth)
        return s.flux * lc
