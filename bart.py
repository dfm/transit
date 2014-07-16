#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["Star", "Planet", "PlanetarySystem"]

import transit
import numpy as np


# Newton's constant in $R_\odot^3 M_\odot^{-1} {days}^{-2}$.
_G = 2945.4625385377644


class Star(object):
    """
    A "star"—in this context—is the massive central body in a
    :class:`PlanetarySystem`.

    :param mass:
        The mass of the star measured in Solar masses. (default: 1.0)

    :param radius:
        The radius of the star measured in Solar radii. (default: 1.0)

    :param flux:
        The un-occulted flux of the star measured in whatever units you feel
        like using. (default: 1.0)

    :param mu1:
        The first quadratic limb darkening coefficient. (default: 0.0)

    :param mu2:
        The second quadratic limb darkening coefficient. (default: 0.0)

    """

    def __init__(self, mass=1.0, radius=1.0, flux=1.0, mu1=0.0, mu2=0.0):
        self.mass = mass
        self.radius = radius
        self.flux = flux
        self.coeffs = (mu1, mu2)

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

    def lnprior(self):
        """
        Compute the log-prior for specified stellar parameterization. This
        function simply returns ``-numpy.inf`` when any parameters are
        un-physical and ``0`` otherwise.

        """
        if self.mass <= 0 or self.radius <= 0:
            return -np.inf
        if not 0 <= self.q1 <= 1 or not 0 <= self.q2 <= 1:
            return -np.inf
        return 0.0


class Planet(object):
    """
    A "planet"—in this context—is a (possibly) massive body orbiting in a
    :class:`PlanetarySystem`. There are several ways to initialize a planet
    and once it has been added to a system using the
    :func:`PlanetarySystem.add_planet` method, they should all be equivalent.
    The orbital elements specified either specify a Keplerian orbit or an
    osculating Keplerian orbit in the context of an N-body system. This object
    includes all sorts of magic for converting between different
    specifications when needed but the base description of the planet and the
    orbit is parameterized by the parameters:

    .. code-block:: python

        (r, mass, a, t0, e, pomega, ix, iy)

    :param r:
        The radius of the planet measured in Solar radii. (default: 0.0)

    :param mass:
        The mass of the planet in Solar masses. (default: 0.0)

    :param a:
        The semi-major axis of the orbit measured in Solar radii. Either this
        parameter or ``period`` must be provided but not both.

    :param period:
        The period of the orbit in days. Either this parameter or ``a`` must
        be provided but not both.

    :param t0:
        The epoch of the orbit in days. (default: 0.0)

    :param e:
        The eccentricity of the orbit. (default: 0.0)

    :param pomega:
        The orientation of the orbital ellipse in radians. (default: 0.0)

    :param ix:
        The relative inclination of the orbital planet along the line-of-sight
        in degrees. This angle is measured differently than you're used to:
        zero degrees is edge on and 90 degrees in face on. This angle will be
        subtracted from the base inclination of the planetary system to get
        the standard measurement of the inclination. Either this parameter
        or ``b`` can be specified but not both. (default: 0.0)

    :param b:
        The mean impact parameter of the orbit measured in stellar radii (not
        Solar radii). Specifically, this impact parameter is defined as

        .. math::
            b = \\frac{a}{R_\star} \\tan i_x

        (default: 0.0)

    :param iy:
        The rotation of the orbital ellipse in the plane of the sky measured
        in radians. Note: this value will not affect the light curve or radial
        velocity values at all. (default: 0.0)

    """

    def __init__(self, r=0.0, mass=0.0, a=None, period=None, t0=0.0, e=0.0,
                 pomega=0.0, ix=None, b=None, iy=0.0):
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
        if not hasattr(self, "planetary_system"):
            raise RuntimeError("You must add this planet to a system "
                               "before getting the period.")

    @property
    def period(self):
        # If we already have a period, return that.
        if self._period is not None:
            return self._period

        # If not, check to make sure that we're already part of a system
        # and then compute the period based on the star's mass.
        self._check_ps()
        mstar = self.planetary_system.star.mass
        a = self._a
        return 2 * np.pi * np.sqrt(a * a * a / _G / (mstar + self.mass))

    @period.setter
    def period(self, P):
        self._check_ps()
        mstar = self.planetary_system.star.mass
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
        rstar = self.planetary_system.star.radius
        return self.a * np.tan(np.radians(self._ix)) / rstar

    @b.setter
    def b(self, b):
        self._check_ps()
        rstar = self.planetary_system.star.radius
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

    def lnprior(self):
        """
        Compute the log-prior for the specified orbital elements. This
        function simply returns ``-numpy.inf`` when any parameters are
        un-physical and ``0`` otherwise.

        """
        if (self.r < 0.0 or self.a <= 0.0 or (not 0.0 <= self.e < 1.0)
           or self.mass < 0.0):
            return -np.inf
        return 0.0


class PlanetarySystem(object):
    """
    A "planetary system" contains a "central" bright :class:`Star` and some
    number (``>= 0``) :class:`Planet` objects orbiting. The default orbits
    are purely Keplerian but sub-classes can include more sophisticated
    solvers.

    :param star:
        A :class:`Star` that specifies the central bright object.

    :para iobs:
        The inclination of the mean orbital plane in degrees. This is
        measured in the standard way with zero inclination meaning face on and
        90 degrees is edge on. (default: 90.0)

    """

    def __init__(self, star, iobs=90.0, hmax=100.0):
        self.star = star
        self.star.planetary_system = self
        self.iobs = iobs
        self.hmax = hmax
        self.planets = []

    @property
    def nplanets(self):
        """
        The number of planets in the system.

        """
        return len(self.planets)

    def add_planet(self, planet):
        """
        Add a :class:`Planet` to the system. This function also sets the
        ``planetary_system`` attribute of the system.

        :param planet:
            The :class:`Planet` to add.

        """
        planet.planetary_system = self
        self.planets.append(planet)

    def _get_pars(self):
        """
        Coerce the planetary parameters into the form needed by the
        :func:`lightcurve` function. This function returns a list of the
        following arrays

        .. code-block:: python

            mass, r, a, t0, e, pomega, ix, iy

        """
        r = [(p.mass, p.r, p.a, p.t0, p.e, p.pomega, p.ix, p.iy)
             for p in self.planets]
        return zip(*r)

    def light_curve(self, t, texp=0.0, tol=0.0, maxdepth=0):
        """
        Get the light curve evaluated at a list of times using the current
        model.

        :param t:
            The times where the light curve should be evaluated (in days).

        """
        if len(self.planets) == 0:
            return self.star.flux + np.zeros_like(t)

        # Coerce the planetary orbital parameters.
        pars = map(np.atleast_1d, self._get_pars())
        mass, r, a, t0, e, pomega, ix, iy = pars

        # Grab some shortcuts to some useful objects.
        s = self.star
        t = np.atleast_1d(t)
        mu1, mu2 = s.coeffs

        # Compute the light curve using the Kepler solver.
        lc = transit.ldlc_kepler(t, mu1, mu2, s.mass, s.radius,
                                 mass, r, a, t0, e, pomega,
                                 np.radians(90.-self.iobs+ix),
                                 np.radians(iy),
                                 texp, tol, maxdepth)
        return s.flux * lc

    def lnprior(self):
        """
        Compute the log prior of the current model. This is simply the sum
        of the log prior of the star and the log priors of the planets.

        """
        lps = self.star.lnprior()
        if not np.isfinite(lps):
            return -np.inf
        lp = np.sum([p.lnprior() for p in self.planets])
        if not np.isfinite(lp):
            return -np.inf
        return lp + lps

    def estimate_t0(self, dataset, grid):
        for ip, planet in enumerate(self.planets):
            max_ll = (-np.inf, 0.0)
            for t0 in grid:
                planet.t0 = t0
                lc = self.light_curve(dataset.time, texp=dataset.texp,
                                      tol=dataset.tol,
                                      maxdepth=dataset.maxdepth)
                ll = dataset.lnlike(lc)
                if ll > max_ll[0]:
                    max_ll = (ll, t0)
            planet.t0 = max_ll[1]
        return max_ll[0]
