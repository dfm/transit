#ifndef _LIMB_DARKENING_H_
#define _LIMB_DARKENING_H_

#include <cmath>

namespace transit {

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
        double lam = 1.0, p2 = p*p;
        if (1+p < z) lam = 0.0;
        else if (z <= 1-p) lam = p2;
        else if (fabs(1-p) < z && z <= 1+p) {
            double z2 = z*z,
                   k1 = acos(0.5*(1-p2+z2)/z),
                   k0 = acos(0.5*(p2+z2-1)/(p*z));
            lam = 1+z2-p2;
            lam = (p2*k0+k1-sqrt(z2-0.25*lam*lam))/M_PI;
        }
        return 1.0 - lam;
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

};

#endif
