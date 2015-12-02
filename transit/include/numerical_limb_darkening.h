#ifndef _NUMERICAL_LIMB_DARKENING_H_
#define _NUMERICAL_LIMB_DARKENING_H_

#include <cmath>
#include "limb_darkening.h"

namespace transit {

class LimbDarkeningLaw {
public:
    virtual ~LimbDarkeningLaw () {};
    virtual double operator () (const double r) const {
        return 1.0;
    };
    virtual double integrate (const double a, const double b) const {
        return M_PI * (b*b - a*a);
    };
    virtual double get_norm () const {
        return M_PI;
    };
};

class QuadraticLaw : public LimbDarkeningLaw {
public:
    QuadraticLaw (const double u1, const double u2) : u1_(u1), u2_(u2) {
        norm_ = integrate(0, 1);
    };
    double operator () (const double r) const {
        double onemmu = 1 - sqrt(1 - r * r);
        return 1 - u1_ * onemmu - u2_ * onemmu * onemmu;
    };
    double integrate (const double a, const double b) const {
        double a2 = a*a,
               b2 = b * b,
               th = 0.5 * 3,
               k1 = 0.5 * (b2 - a2) * (1 - u1_ - 2 * u2_),
               k2 = (u1_ + 2 * u2_) * (pow(1 - a2, th) - pow(1 - b2, th))/3.0,
               k3 = 0.25 * u2_ * (b2 * b2 - a2 * a2);
        return 2 * M_PI * (k1 + k2 + k3);
    };
    double get_norm () const { return norm_; };

private:
    double u1_, u2_, norm_;
};

class NumericalLimbDarkening : public LimbDarkening {
public:
    NumericalLimbDarkening (LimbDarkeningLaw* law, double tol,
                            int maxiter, int iterstep)
        : tol_(tol), maxiter_(maxiter), iterstep_(iterstep), law_(law) {};

    template <typename T>
    T operator () (const T* const params, const T& p, const T& z) const {
        T q1 = params[0], q2 = 2.0 * params[1],
          u1 = q1 * q2, u2 = q1 * (1.0 - q2);

        z = fabs(z);
        double zmp = fmin(fmax(z - p, 0.0), 1.0),
               zpp = fmin(fmax(z + p, 0.0), 1.0);
        if (zpp <= zmp) return T(1.0);

        // Iterate the integration to convergence.
        T norm = law_->get_norm(),
          value = integrate(zmp, zpp, p, z, 1) / norm,
          value0 = value;
        for (int i = 1; i <= maxiter_; i += iterstep_) {
            value = integrate(zmp, zpp, p, z, i+1) / norm;
            if (value != value) {
                return 1.0 - value0;
            }
            if (fabs(value - value0) < tol_) {
                return 1.0 - value;
            }
            value0 = value;
        }
        std::cerr << "No convergence\n";
        return 1.0 - value;
    };

    template <typename T>
    T integrate (const T& a, const T& b, const T& p, const T& z, int N) const {
        T result = T(0.0), r0 = a, oa0 = T(0.0), area0 = a*a,
          r, dr = (b - a) / N, oa, area;
        for (int i = 1; i <= N; ++i) {
            // Compute the areas at this radius.
            r = a + i*dr;
            area = r*r;
            oa = occ_area(r, p, z);
            if (oa != oa)
                std::cout << "DUDEFACE " << r << " " << p << " " << z << " " << oa << std::endl;

            // Update the result.
            result += (oa - oa0) * law_->integrate(r0, r) / (area - area0);
            if (result != result) {
                std::cout << r << " " << p << " " << z << " " << oa << " " << oa0 << " blurghleshorts\n";
            }

            // Cache the results of this step.
            r0 = r;
            area0 = area;
            oa0 = oa;
        }
        return result;
    };

private:
    int miniter_, maxiter_, iterstep_;
    double tol_;
    LimbDarkeningLaw* law_;
};

};

#endif
