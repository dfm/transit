#ifndef _LIMB_DARKENING_H_
#define _LIMB_DARKENING_H_

#include <cmath>

namespace transit {

inline double _acos (const double& x) {
    return std::acos(fmax(fmin(x, 1.0), -1.0));
}

inline double occ_area (const double& r0, const double& p, const double& b)
{
    if (b > r0 + p) return 0.0;
    else if (fabs(r0-p) < b && b <= r0+p) {
        double r2 = r0*r0, p2 = p*p, b2 = b*b, k1, k2, k3;
        k1 = _acos(0.5 * (b2 + p2 - r2) / b / p);
        k2 = _acos(0.5 * (b2 + r2 - p2) / b / r0);
        k3 = sqrt((p+r0-b) * (b+p-r0) * (b-p+r0) * (b+r0+p));
        return (p2 * k1 + r2 * k2 - 0.5 * k3) / M_PI;
    } else if (b <= r0 - p) return p * p;
    return r0 * r0;
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
        return 1.0 - occ_area(1.0, p, z);
    };
};

class QuadraticLimbDarkening : public LimbDarkening {
public:
    QuadraticLimbDarkening () {};
    ~QuadraticLimbDarkening () {};
    QuadraticLimbDarkening (double u1, double u2)
        : u1_(u1), u2_(u2) {};
    double operator () (double p, double z) const;

private:
    double u1_, u2_;
};

};

#endif
