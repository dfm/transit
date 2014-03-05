#ifndef _LIMB_DARKENING_H_
#define _LIMB_DARKENING_H_

#include <cmath>
#include "autodiff.h"

namespace transit {

template <typename T>
T occ_area (const double& r0, const T& p, const T& b)
{
    if (b >= r0 + p) return T(0.0);
    else if (b <= r0 - p) return T(p * p);
    else if (b <= p - r0) return T(r0 * r0);

    double r2 = r0*r0;
    T p2 = p*p, b2 = b*b, k1, k2, k3;
    k1 = acos(T(0.5) * (b2 + p2 - r2) / b / p);
    k2 = acos(T(0.5) * (b2 + r2 - p2) / b / r0);
    k3 = sqrt((p+r0-b) * (b+p-r0) * (b-p+r0) * (b+r0+p));
    return (p2 * k1 + r2 * k2 - T(0.5) * k3) / T(M_PI);
}

class LimbDarkening {
public:
    virtual ~LimbDarkening () {};
    virtual double operator () (double p, double z) const {
        return 1.0;
    };
};

class GeometricLimbDarkening : public LimbDarkening {
public:
    double operator () (double p, double z) const {
        return 1.0 - occ_area<double>(1.0, p, z);
    };
};

class QuadraticLimbDarkening : public LimbDarkening {
public:
    QuadraticLimbDarkening (double u1, double u2)
        : u1_(u1), u2_(u2) {};
    double operator () (double p, double z) const;

private:
    double u1_, u2_;
};

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
    double operator () (double p, double z) const {
        z = fabs(z);
        double zmp = fmin(fmax(z - p, 0), 1),
               zpp = fmin(fmax(z + p, 0), 1);
        if (zpp <= zmp) return 1.0;

        // Iterate the integration to convergence.
        double norm = law_->get_norm(),
               value = integrate(zmp, zpp, p, z, 1) / norm,
               value0 = value;
        for (int i = 1; i <= maxiter_; i += iterstep_) {
            value = integrate(zmp, zpp, p, z, i+1) / norm;
            if (fabs(value - value0) < tol_) {
                std::cout << "Converged after " << i << " iterations\n";
                return 1.0 - value;
            }
            value0 = value;
        }
        std::cout << "No convergence\n";
        return 1.0 - value;
    };

    double integrate (double a, double b, double p, double z, int N) const {
        double result = 0.0, r0 = a, oa0 = 0.0, area0 = a*a,
               r, dr = (b - a) / N, oa, area;
        for (int i = 1; i <= N; ++i) {
            // Compute the areas at this radius.
            r = a + i*dr;
            area = r*r;
            oa = occ_area<double>(r, p, z);

            // Update the result.
            result += (oa - oa0) * law_->integrate(r0, r) / (area - area0);

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

class NumericalLimbDarkening2 : public LimbDarkening {
public:
    NumericalLimbDarkening2 (LimbDarkeningLaw* law, double tol,
                             int maxiter, int iterstep)
        : tol_(tol), maxiter_(maxiter), iterstep_(iterstep), law_(law) {};
    double operator () (double p, double z) const {
        z = fabs(z);
        double zmp = fmin(fmax(z - p, 0), 1),
               zpp = fmin(fmax(z + p, 0), 1);
        if (zpp <= zmp) return 1.0;

        return integrate(0, 1, p, z, 10000);

        // // Iterate the integration to convergence.
        // double norm = law_->get_norm(),
        //        value = integrate(zmp, zpp, p, z, 1) / norm,
        //        value0 = value;
        // for (int i = 1; i <= maxiter_; i += iterstep_) {
        //     value = integrate(zmp, zpp, p, z, i+1) / norm;
        //     if (fabs(value - value0) < tol_) {
        //         std::cout << "Converged after " << i << " iterations\n";
        //         return 1.0 - value;
        //     }
        //     value0 = value;
        // }
        // std::cout << "No convergence\n";
        // return 1.0 - value;
    };

    double integrate (double a, double b, double p, double z, int N) const {
        double dr = (b - a) / N, result = 0.0, norm = 0.0;
        for (int i = 0; i < N; ++i) {
            AutoDiff<1> r(a + (i+0.5)*dr, 0),
                        f = r * r * (1 - occ_area(1.0, p/r, z/r));
            double Ir = (*law_)(r.a);
            result += f.v(0) * Ir;
            norm += 2 * r.a * Ir;
            std::cout << Ir << " " << r.a << " " << f.v(0) << std::endl;
        }
        std::cout << result * dr << " " << norm * dr << std::endl;
        return result / norm;
    };

private:
    int miniter_, maxiter_, iterstep_;
    double tol_;
    LimbDarkeningLaw* law_;
};

};

#endif
