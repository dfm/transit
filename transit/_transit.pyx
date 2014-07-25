# distutils: language = c++

cdef extern from "limb_darkening.h" namespace "transit":

    cdef cppclass LimbDarkening:
        pass

    cdef cppclass QuadraticLimbDarkening:
        QuadraticLimbDarkening () except +
        QuadraticLimbDarkening (double, double) except +


cdef extern from "kepler.h" namespace "transit":

    cdef cppclass KeplerSolver[L]:
        KeplerSolver (L, double, double) except +
        void position (const double t, const int i, double pos[])
        int add_body (LimbDarkening*, double occ, double m, double r,
                      double a, double t0, double e, double pomega, double ix,
                      double iy)


cdef class PythonKeplerSolver:
    cdef KeplerSolver[QuadraticLimbDarkening] *thisptr

    def __cinit__(self, double u1, double u2, double mstar, double rstar):
        cdef QuadraticLimbDarkening ld = QuadraticLimbDarkening(u1, u2)
        self.thisptr = new KeplerSolver[QuadraticLimbDarkening](ld, mstar,
                                                                rstar)

    def __dealloc__(self):
        del self.thisptr
