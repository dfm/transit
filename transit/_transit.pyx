# distutils: language = c++
from __future__ import division

cimport cython
from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

DTYPE_u = np.uint32
ctypedef np.uint32_t DTYPE_u_t


cdef extern from "quad.h" namespace "transit":
    cdef cppclass QuadraticLimbDarkening:
        QuadraticLimbDarkening () except +


cdef extern from *:
    ctypedef int simple_type "7"
    ctypedef int kepler_0 "5"
    ctypedef int kepler_1 "5+7"
    ctypedef int kepler_2 "5+7*2"
    ctypedef int kepler_3 "5+7*3"
    ctypedef int kepler_4 "5+7*4"
    ctypedef int kepler_5 "5+7*5"

    ctypedef int ttvfaster_2 "1+7*2"
    ctypedef int ttvfaster_3 "1+7*3"
    ctypedef int ttvfaster_4 "1+7*4"
    ctypedef int ttvfaster_5 "1+7*5"


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
        void* reparameterize (const void* const params)


cdef extern from "simple.h" namespace "transit":
    cdef cppclass SimpleSolver[L]:
        SimpleSolver ()
        int get_status () const
        void* reparameterize (const void* const params)

cdef extern from "ttvfaster.h" namespace "transit::ttvfaster":
    unsigned compute_ntransits (
        unsigned n_planets,
        double* params,
        double t0,
        double tf,
        unsigned* n_transits,
        unsigned* starts)
    void compute_times (
        unsigned n_planets,
        double* params,
        double t0,
        unsigned m_max,
        unsigned* n_transits,
        unsigned* starts,
        double* times)
    void compute_grad_times[N] (
        unsigned n_planets,
        double* params,
        double t0,
        unsigned m_max,
        unsigned* n_transits,
        unsigned* starts,
        double* times,
        double* grad_times)


cdef class CythonSolver:

    def simple_light_curve(self,
                           np.ndarray[DTYPE_t] params,
                           np.ndarray[DTYPE_t] t,
                           double texp, double tol, int maxdepth):
        cdef double* reparams
        cdef Integrator[SimpleSolver[QuadraticLimbDarkening] ] integrator = \
            Integrator[SimpleSolver[QuadraticLimbDarkening] ] (tol, maxdepth)

        # Compute the light curve.
        cdef np.ndarray[DTYPE_t] lc = np.empty(t.shape[0], dtype=DTYPE)
        cdef int flag = integrator.integrate(<double*>params.data,
                                             t.shape[0],
                                             <double*>t.data,
                                             texp,
                                             <double*>lc.data)

        return lc

    def simple_gradient(self,
                        np.ndarray[DTYPE_t] params,
                        np.ndarray[DTYPE_t] t,
                        double texp, double tol, int maxdepth):
        cdef double* reparams
        cdef Integrator[SimpleSolver[QuadraticLimbDarkening] ] integrator = \
            Integrator[SimpleSolver[QuadraticLimbDarkening] ] (tol, maxdepth)

        # Compute the light curve.
        cdef np.ndarray[DTYPE_t] lc = np.empty(t.shape[0], dtype=DTYPE)
        cdef np.ndarray[DTYPE_t, ndim=2] gradient = \
            np.empty((t.shape[0], params.shape[0]), dtype=DTYPE)

        cdef int flag

        flag = integrator.gradient[simple_type](
            <double*>params.data, t.shape[0], <double*>t.data, texp,
            <double*>lc.data, <double*>gradient.data)

        if flag:
            raise RuntimeError("the solver failed with code={0}".format(flag))

        return lc, gradient

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

    def approx_times(self,
                     unsigned n_planets,
                     np.ndarray[DTYPE_t] params,
                     double tmin,
                     double tmax,
                     unsigned m_max):
        cdef unsigned i, j, n, ind
        cdef np.ndarray[DTYPE_u_t] n_transits = np.empty(n_planets, dtype=DTYPE_u)
        cdef np.ndarray[DTYPE_u_t] starts = np.empty(n_planets, dtype=DTYPE_u)

        cdef unsigned ntot = compute_ntransits (
            n_planets, <double*>params.data, tmin, tmax,
            <unsigned*>n_transits.data, <unsigned*>starts.data)

        cdef np.ndarray[DTYPE_t] times = np.empty(ntot, dtype=DTYPE)
        compute_times (
            n_planets, <double*>params.data, tmin, m_max,
            <unsigned*>n_transits.data, <unsigned*>starts.data,
            <double*>times.data
        )

        del starts
        return n_transits, times

    def grad_approx_times(self,
                          unsigned n_planets,
                          np.ndarray[DTYPE_t] params,
                          double tmin,
                          double tmax,
                          unsigned m_max):
        cdef unsigned i, j, n, ind
        cdef np.ndarray[DTYPE_u_t] n_transits = np.empty(n_planets, dtype=DTYPE_u)
        cdef np.ndarray[DTYPE_u_t] starts = np.empty(n_planets, dtype=DTYPE_u)

        cdef unsigned ntot = compute_ntransits (
            n_planets, <double*>params.data, tmin, tmax,
            <unsigned*>n_transits.data, <unsigned*>starts.data)

        cdef np.ndarray[DTYPE_t] times = np.empty(ntot, dtype=DTYPE)
        cdef np.ndarray[DTYPE_t, ndim=2] grad_times = \
            np.empty((ntot, 1 + n_planets*7), dtype=DTYPE)

        if n_planets == 2:
            compute_grad_times[ttvfaster_2] (
                n_planets, <double*>params.data, tmin, m_max,
                <unsigned*>n_transits.data, <unsigned*>starts.data,
                <double*>times.data, <double*>grad_times.data
            )
        else:
            assert False

        del starts
        return n_transits, times, grad_times
