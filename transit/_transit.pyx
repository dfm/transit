# distutils: language = c++
from __future__ import division

cimport cython
from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


cdef extern from "quad.h" namespace "transit":
    cdef cppclass QuadraticLimbDarkening:
        QuadraticLimbDarkening () except +


cdef extern from *:
    ctypedef int kepler_0 "5"
    ctypedef int kepler_1 "5+8"
    ctypedef int kepler_2 "5+8*2"
    ctypedef int kepler_3 "5+8*3"
    ctypedef int kepler_4 "5+8*4"
    ctypedef int kepler_5 "5+8*5"
    ctypedef int kepler_6 "5+8*6"
    ctypedef int kepler_7 "5+8*7"
    ctypedef int kepler_8 "5+8*8"
    ctypedef int kepler_9 "5+8*9"


cdef extern from "integrator.h" namespace "transit":
    cdef cppclass Integrator[S]:
        Integrator ()
        Integrator (double tol, int maxdepth)
        int get_status () const
        S* get_solver ()
        int integrate (const double* params, int nt, const double* t,
                       double texp, double* lc)
        int gradient[N](const double* params, int nt, const double* t,
                        double texp, double* lc, double* gradient)


cdef extern from "kepler.h" namespace "transit":
    cdef cppclass KeplerSolver[L]:
        KeplerSolver ()
        int get_status () const
        void set_n_body (int n)
        int get_n_body () const
        void position (const double t, const int i, double pos[])
        void* reparameterize (const void* const params)


# cdef extern from "solver.h" namespace "transit":

#     cdef cppclass SimpleSolver[L]:
#         SimpleSolver ()
#         SimpleSolver (L* ld, double period, double t0, double duration,
#                       double ror, double impact) except +


# cdef class PythonSimpleSolver:
#     cdef SimpleSolver[QuadraticLimbDarkening] *thisptr

#     def __cinit__(self, double u1, double u2, double period, double t0,
#                   double duration, double ror, double impact):
#         cdef QuadraticLimbDarkening* ld = new QuadraticLimbDarkening(u1, u2)
#         self.thisptr = new SimpleSolver[QuadraticLimbDarkening](
#             ld, period, t0, duration, ror, impact)

#     def __dealloc__(self):
#         del self.thisptr

#     def light_curve(self, double f0, np.ndarray[DTYPE_t] t, double texp,
#                     double tol, int maxdepth):
#         cdef int i, N = t.shape[0]

#         # Define the integrator.
#         cdef Integrator[SimpleSolver[QuadraticLimbDarkening] ] integrator = \
#             Integrator[SimpleSolver[QuadraticLimbDarkening] ] \
#             (self.thisptr, tol, maxdepth)

#         # Compute the integrated light curve.
#         cdef np.ndarray[DTYPE_t] lam = np.zeros(N, dtype=DTYPE)
#         for i in range(N):
#             lam[i] = f0 * integrator.integrate(t[i], texp)
#             if integrator.get_status():
#                 raise RuntimeError("Integrator failed with status {0}"
#                                    .format(integrator.get_status()))

#         return lam


cdef class CythonSolver:

    def kepler_light_curve(self,
                           int n_body,
                           np.ndarray[DTYPE_t] params,
                           np.ndarray[DTYPE_t] t,
                           double texp, double tol, int maxdepth):
        cdef double* reparams
        cdef Integrator[KeplerSolver[QuadraticLimbDarkening] ] integrator = \
            Integrator[KeplerSolver[QuadraticLimbDarkening] ] (tol, maxdepth)
        integrator.get_solver().set_n_body(n_body)

        # Compute the light curve.
        cdef np.ndarray[DTYPE_t] lc = np.empty(t.shape[0], dtype=DTYPE)
        cdef int flag = integrator.integrate(<double*>params.data,
                                             t.shape[0],
                                             <double*>t.data,
                                             texp,
                                             <double*>lc.data)

        return lc

    def kepler_gradient(self,
                        int n_body,
                        np.ndarray[DTYPE_t] params,
                        np.ndarray[DTYPE_t] t,
                        double texp, double tol, int maxdepth):
        cdef double* reparams
        cdef Integrator[KeplerSolver[QuadraticLimbDarkening] ] integrator = \
            Integrator[KeplerSolver[QuadraticLimbDarkening] ] (tol, maxdepth)
        integrator.get_solver().set_n_body(n_body)

        # Compute the light curve.
        cdef np.ndarray[DTYPE_t] lc = np.empty(t.shape[0], dtype=DTYPE)
        cdef np.ndarray[DTYPE_t, ndim=2] gradient = \
            np.empty((t.shape[0], params.shape[0]), dtype=DTYPE)

        cdef int flag

        if n_body == 0:
            flag = integrator.gradient[kepler_0](
                <double*>params.data, t.shape[0], <double*>t.data, texp,
                <double*>lc.data, <double*>gradient.data)
        elif n_body == 1:
            flag = integrator.gradient[kepler_1](
                <double*>params.data, t.shape[0], <double*>t.data, texp,
                <double*>lc.data, <double*>gradient.data)
        elif n_body == 2:
            flag = integrator.gradient[kepler_2](
                <double*>params.data, t.shape[0], <double*>t.data, texp,
                <double*>lc.data, <double*>gradient.data)
        elif n_body == 3:
            flag = integrator.gradient[kepler_3](
                <double*>params.data, t.shape[0], <double*>t.data, texp,
                <double*>lc.data, <double*>gradient.data)
        elif n_body == 4:
            flag = integrator.gradient[kepler_4](
                <double*>params.data, t.shape[0], <double*>t.data, texp,
                <double*>lc.data, <double*>gradient.data)
        elif n_body == 5:
            flag = integrator.gradient[kepler_5](
                <double*>params.data, t.shape[0], <double*>t.data, texp,
                <double*>lc.data, <double*>gradient.data)
        else:
            raise ValueError("You're a maniac. Don't take the gradient of a "
                             "Kepler light curve with more than 5 planets!")

        if flag:
            raise RuntimeError("the solver failed with code={0}".format(flag))

        return lc, gradient
