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
        double radial_velocity (const double t, const int i)


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
        cdef int info
        info = self.thisptr.add_body (ld, occ, m, r, a, t0, e, pomega, ix, iy)
        if info:
            raise RuntimeError("Couldn't add body.")

    @cython.boundscheck(False)
    def position(self, np.ndarray[DTYPE_t, ndim=1] t):
        cdef unsigned int i, j, k
        cdef unsigned int n = t.shape[0], K = self.thisptr.nbodies()

        cdef double value[3]
        cdef np.ndarray[DTYPE_t, ndim=3] pos = np.zeros([n, K, 3],
                                                        dtype=DTYPE)
        for i in range(n):
            for j in range(K):
                self.thisptr.position(t[i], j, value)
                if self.status:
                    raise RuntimeError("Kepler solver failed with status {0}"
                                    .format(self.status))
                for k in range(3):
                    pos[i, j, k] = value[k]
        return pos

    @cython.boundscheck(False)
    def velocity(self, np.ndarray[DTYPE_t, ndim=1] t):
        cdef unsigned int i, j, k
        cdef unsigned int n = t.shape[0], K = self.thisptr.nbodies()

        cdef double value[3]
        cdef np.ndarray[DTYPE_t, ndim=3] vel = np.zeros([n, K, 3],
                                                        dtype=DTYPE)
        for i in range(n):
            for j in range(K):
                self.thisptr.velocity(t[i], j, value)
                if self.status:
                    raise RuntimeError("Kepler solver failed with status {0}"
                                    .format(self.status))
                for k in range(3):
                    vel[i, j, k] = value[k]
        return vel

    def light_curve(self, double f0, np.ndarray[DTYPE_t] t, double texp,
                    double tol, int maxdepth):
        cdef int i, N = t.shape[0]

        # Define the integrator.
        cdef Integrator[KeplerSolver[QuadraticLimbDarkening] ] integrator = \
                Integrator[KeplerSolver[QuadraticLimbDarkening] ] \
                    (self.thisptr, tol, maxdepth)

        # Compute the integrated light curve.
        cdef np.ndarray[DTYPE_t] lam = np.zeros(N, dtype=DTYPE)
        for i in range(N):
            lam[i] = f0 * integrator.integrate(t[i], texp)
            if integrator.get_status():
                raise RuntimeError("Integrator failed with status {0}"
                                   .format(integrator.get_status()))

        return lam
