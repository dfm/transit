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

    cdef cppclass QuadraticLimbDarkening:
        QuadraticLimbDarkening () except +
        QuadraticLimbDarkening (double, double) except +


cdef extern from "kepler.h" namespace "transit":

    cdef cppclass KeplerSolver[L]:
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
            for j in range(3):
                vel[i, j] = value[j]
        return vel
