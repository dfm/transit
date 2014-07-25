#ifndef _KEPLER_H_
#define _KEPLER_H_

#include <cfloat>
#include <vector>
#include <cstddef>

#include "autodiff.h"
#include "limb_darkening.h"

#define KEPLER_MAX_ITER 200
#define KEPLER_CONV_TOL 1.48e-10

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//
// Newton's constant in $R_\odot^3 M_\odot^{-1} {days}^{-2}$.
//
#define G_GRAV 2945.4625385377644

namespace transit {

double ecc_to_mean_anomaly (double psi, double e) {
    return psi - e * sin(psi);
}

double mean_to_ecc_anomaly (double M, double e, int *info) {
    int i;
    double wt0, psi0, psi = 0.0;

    // Check for un-physical parameters.
    if (e < 0 || e >= 1) {
        *info = 2;
        return 0.0;
    }

    *info = 0;
    wt0 = fmod(M, 2 * M_PI);
    psi0 = wt0;

    for (i = 0; i < KEPLER_MAX_ITER; ++i) {
        // Compute the function value and the first two derivatives.
        double spsi0 = sin(psi0),
               f = psi0 - e * spsi0 - wt0,
               fp = 1.0 - e * cos(psi0),
               fpp = e * spsi0;

        // Take a second order step.
        psi = psi0 - 2.0 * f * fp / (2.0 * fp * fp - f * fpp);

        // Accept?
        if (fabs(psi - psi0) <= KEPLER_CONV_TOL) return psi;

        // Save as the previous step.
        psi0 = psi;
    }

    // If we get here, we didn't ever converge.
    *info = 1;
    return psi;
}

Jet<double> mean_to_ecc_anomaly (Jet<double> M, double e, int *info) {
    double psi = mean_to_ecc_anomaly(M.a, e, info);
    return Jet<double>(psi, M.v / (1 - e * cos(psi)));
}

double ecc_to_true_anomaly (double E, double e, int* info)
{
    double denom = 1. - e * cos(E);
    if (denom <= 0.0) {
        *info = -1;
        return E;
    }
    double f = acos((cos(E) - e) / denom);
    return f * sgn(sin(f)) * sgn(sin(E));
}

template <typename T>
int solve_kepler (const T& manom, T* pos,
                  const double a, const double e,
                  const double sp, const double cp,
                  const double six, const double cix,
                  const double siy, const double ciy)
{
    int info;
    T psi = mean_to_ecc_anomaly (manom, e, &info);

    // Did solve fail?
    if (info) return info;

    // Compute the true anomaly.
    T cpsi = cos(psi), spsi = sin(psi), d = 1.0 - e * cpsi;
    if (d == 0.0) return 3;
    T cth = (cpsi - e) / d, sth = sqrt(1.0 - cth * cth);

    // Compute the radius and coordinates in the plane.
    d = (1.0 + e * cth);
    if (d == 0.0) return 3;
    T r = a * (1.0 - e * e) / d,
      x = r * cth,
      y = r * sth * T(sgn(spsi));

    // Rotate by pomega.
    T xp =  x * cp + y * sp,
      yp = -x * sp + y * cp;

    // Compute the positions.
    pos[0] = xp * cix;
    pos[1] = yp * ciy - xp * six * siy;
    pos[2] = yp * siy + xp * six * ciy;
    return 0;
}

template <class L>
class KeplerSolver {

public:

    KeplerSolver () {};
    KeplerSolver (L* ld, double mstar, double rstar)
        : ld_(ld), mstar_(mstar), rstar_(rstar)
    {
        status_ = 0;
        f0_ = 1.0;
    };

    ~KeplerSolver () {
        delete ld_;
    };

    int get_status () const { return status_; };
    int nbodies () const { return mp_.size(); };

    int add_body (LimbDarkening* ld, double occ, double m, double r, double a,
                  double t0, double e, double pomega, double ix, double iy)
    {
        if (e < 0.0 || fabs(e - 1.0) <= DBL_EPSILON) return 1;

        pld_.push_back(ld);
        occ_.push_back(occ);
        mp_.push_back(m);
        r_.push_back(r);
        ror_.push_back(r/rstar_);
        iror_.push_back(rstar_/r);
        a_.push_back(a);
        t0_.push_back(t0);
        e_.push_back(e);
        pomega_.push_back(pomega);
        ix_.push_back(ix);
        iy_.push_back(iy);

        // Pre-compute some constant factors.
        double period = 2 * M_PI * sqrt(a*a*a/G_GRAV/(mstar_+m)),
               psi0 = 2 * atan2(tan(0.5 * pomega), sqrt((1 + e) / (1 - e)));
        periods_.push_back(period);
        dmanomdt_.push_back(2 * M_PI / period);
        t1s_.push_back(t0-0.5*period*(psi0 - e * sin(psi0)) / M_PI);
        cpom_.push_back(cos(pomega));
        spom_.push_back(sin(pomega));
        cix_.push_back(cos(ix));
        six_.push_back(sin(ix));
        ciy_.push_back(cos(iy));
        siy_.push_back(sin(iy));

        // Radial velocity term.
        K_.push_back(pow(2* M_PI * G_GRAV / period, 1./3) * m * sin(ix)
                     / (pow(m + mstar_, 2 / 3.) * sqrt(1 - e*e)));

        return 0;
    };

    void position (const double t, const int i, double pos[]) {
        double manom = dmanomdt_[i] * (t - t1s_[i]);
        status_ = solve_kepler (manom, pos, a_[i], e_[i], spom_[i], cpom_[i],
                                six_[i], cix_[i], siy_[i], ciy_[i]);
    };

    void velocity (const double t, int i, double vel[]) {
        Jet<double> manom(dmanomdt_[i] * (t - t1s_[i]), dmanomdt_[i]),
                    res[3];

        status_ = solve_kepler (manom, res, a_[i], e_[i], spom_[i], cpom_[i],
                                six_[i], cix_[i], siy_[i], ciy_[i]);

        vel[0] = res[0].v;
        vel[1] = res[1].v;
        vel[2] = res[2].v;
    };

    virtual double operator () (const double t) {
        int i, l = mp_.size();
        double z, lam = 1.0, lp = 0.0, pos[3] = {0, 0, 0};
        for (i = 0; i < l; ++i) {
            // Solve Kepler's equation for the position of the body.
            position(t, i, pos);

            // Find the impact parameter in the plane of the sky.
            z = sqrt(pos[1]*pos[1] + pos[2]*pos[2]);

            // Is the body transiting or occulting.
            if (pos[0] > 0.0) {
                // If transiting, use the stellar limb darkening profile.
                lam *= (*ld_)(ror_[i], z/rstar_);
                lp += occ_[i];
            } else if (occ_[i] > 0.0 && pld_[i] != NULL) {
                // If occulting, use the planet's limb darkening profile.
                double l = (*(pld_[i]))(iror_[i], z/r_[i]);
                lp += occ_[i] * l;
            }
        }
        return lam + lp;
    };

    // double radial_velocity (const double t, int i) {
    //     int info;
    //     double e = e_[i], pomega = pomega_[i], spomega = spom_[i],
    //            dmanomdt = dmanomdt_[i], t1 = t1s_[i],
    //            manom = dmanomdt * (t - t1),
    //            psi = mean_to_ecc_anomaly (manom, e, &info),
    //            f = ecc_to_true_anomaly (psi, e),
    //            K = K_[i];
    //     return K * (sin(f+pomega) + e * spomega);
    // }


protected:

    int status_;
    L* ld_;
    double mstar_, rstar_, f0_;
    std::vector<LimbDarkening*> pld_;
    std::vector<double> occ_, mp_, r_, ror_, iror_, a_, t0_, e_, pomega_, ix_,
                        iy_, periods_, dmanomdt_, t1s_, cpom_, spom_,
                        cix_, six_, ciy_, siy_, K_;

};

};

#endif
