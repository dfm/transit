# distutils: language = c++
from __future__ import division

cimport cython

import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

cdef extern from "limb_darkening.h" namespace "transit":

    cdef cppclass LimbDarkening:
        pass

    cdef cppclass GeometricLimbDarkening(LimbDarkening):
        pass

    cdef cppclass QuadraticLimbDarkening(LimbDarkening):
        QuadraticLimbDarkening () except +
        QuadraticLimbDarkening (double, double) except +


cdef extern from "integrator.h" namespace "transit":

    cdef cppclass Integrator[S]:
        Integrator ()
        Integrator (S*, double, int)
        int get_status () const
        double integrate (double f0, double t, double texp, int depth)
        double integrate (double t, double texp)


cdef extern from "kepler.h" namespace "transit":

    cdef cppclass KeplerSolver[L]:
        KeplerSolver ()
        KeplerSolver (L*, double, double) except +
        int get_status () const
        int nbodies () const
        int add_body (LimbDarkening*, double occ, double m, double r,
                      double a, double t0, double e, double pomega, double ix,
                      double iy)
        void position (const double t, const int i, double pos[])
        void velocity (const double t, const int i, double vel[])


cdef class PythonKeplerSolver:
    cdef KeplerSolver[QuadraticLimbDarkening] *thisptr

    def __cinit__(self, double u1, double u2, double mstar, double rstar):
        cdef QuadraticLimbDarkening* ld = new QuadraticLimbDarkening(u1, u2)
        self.thisptr = new KeplerSolver[QuadraticLimbDarkening](ld, mstar,
                                                                rstar)

    def __dealloc__(self):
        del self.thisptr

    property status:
        def __get__(self): return self.thisptr.get_status()

    def __len__(self):
        return self.thisptr.nbodies()

    def add_body(self, double occ, double m, double r, double a, double t0,
                 double e, double pomega, double ix, double iy):
        cdef LimbDarkening* ld = new LimbDarkening()
        self.thisptr.add_body (ld, occ, m, r, a, t0, e, pomega, ix, iy)

    @cython.boundscheck(False)
    def get_position(self, np.ndarray[DTYPE_t, ndim=1] t, int ind):
        cdef unsigned int i, j
        cdef double value[3]
        cdef int n = t.shape[0]
        cdef np.ndarray[DTYPE_t, ndim=2] pos = np.zeros([n, 3], dtype=DTYPE)
        for i in range(n):
            self.thisptr.position(t[i], ind, value)
            if self.status:
                raise RuntimeError("Kepler solver failed with status {0}"
                                   .format(self.status))
            for j in range(3):
                pos[i, j] = value[j]
        return pos

    @cython.boundscheck(False)
    def get_velocity(self, np.ndarray[DTYPE_t, ndim=1] t, int ind):
        cdef unsigned int i, j
        cdef double value[3]
        cdef int n = t.shape[0]
        cdef np.ndarray[DTYPE_t, ndim=2] vel = np.zeros([n, 3], dtype=DTYPE)
        for i in range(n):
            self.thisptr.velocity(t[i], ind, value)
            if self.status:
                raise RuntimeError("Kepler solver failed with status {0}"
                                   .format(self.status))
            for j in range(3):
                vel[i, j] = value[j]
        return vel


def ldlc_kepler (np.ndarray[DTYPE_t] t,
                 double u1, double u2, double mstar, double rstar,
                 np.ndarray[DTYPE_t] occ, np.ndarray[DTYPE_t] m,
                 np.ndarray[DTYPE_t] r, np.ndarray[DTYPE_t] a,
                 np.ndarray[DTYPE_t] t0, np.ndarray[DTYPE_t] e,
                 np.ndarray[DTYPE_t] pomega, np.ndarray[DTYPE_t] ix,
                 np.ndarray[DTYPE_t] iy,
                 double texp, double tol, int maxdepth):
    cdef int i, info, nplanets = occ.shape[0], N = t.shape[0]

    cdef GeometricLimbDarkening pld
    cdef QuadraticLimbDarkening* ld = new QuadraticLimbDarkening(u1, u2)
    cdef KeplerSolver[QuadraticLimbDarkening]* solver = \
            new KeplerSolver[QuadraticLimbDarkening] (ld, mstar, rstar)

    # Add the planets to the system.
    for i in range(nplanets):
        info = solver.add_body(&pld, occ[i], m[i], r[i], a[i], t0[i], e[i],
                               pomega[i], ix[i], iy[i])
        if info:
            del solver
            raise RuntimeError("Couldn't add body {0}".format(i))

    cdef Integrator[KeplerSolver[QuadraticLimbDarkening] ] integrator = \
            Integrator[KeplerSolver[QuadraticLimbDarkening] ](solver, tol,
                                                              maxdepth)

    # Compute the integrated light curve.
    cdef np.ndarray[DTYPE_t] lam = np.zeros(N, dtype=DTYPE)
    for i in range(N):
        lam[i] = integrator.integrate(t[i], texp)
        if integrator.get_status():
            del solver
            raise RuntimeError("Integrator failed with status {0}"
                               .format(integrator.get_status()))

    return lam
