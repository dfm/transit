#ifndef _KEPLER_H_
#define _KEPLER_H_

#include <cmath>
#include <cfloat>
#include <vector>

#include "limb_darkening.h"

#define KEPLER_MAX_ITER 200
#define KEPLER_CONV_TOL 1.48e-10

//
// Newton's constant in $R_\odot^3 M_\odot^{-1} {days}^{-2}$.
//
#define G_GRAV 2945.4625385377644

namespace transit {

double kepler_ecc_to_mean_anomaly (double psi, double e) {
    return psi - e * sin(psi);
}

double kepler_mean_to_ecc_anomaly (double M, double e, int *info) {
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

template <class L>
class KeplerSolver {

public:

    KeplerSolver (L ld, double mstar, double rstar)
        : ld_(ld), mstar_(mstar), rstar_(rstar)
    {
        status_ = 0;
    };

    int get_status () const { return status_; };

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

        // // Doppler beaming, etc.
        // zp_.push_back(zp);
        // ae_.push_back(ae);
        // lag_.push_back(lag);
        // ab_.push_back(ab);
        // ar_.push_back(ar);

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

        return 0;
    };

    void solve (double t, int i, double pos[]) {
        int info = 0, sgn;
        double e = e_[i], a = a_[i],
               dmanomdt = dmanomdt_[i],
               t1 = t1s_[i],
               cpomega = cpom_[i], spomega = spom_[i],
               cix = cix_[i], six = six_[i], ciy = ciy_[i], siy = siy_[i],
               manom, psi, cpsi, spsi, d, cth, sth, r, x, y, xp, yp, xsx,
               denom;

        manom = dmanomdt * (t - t1);
        psi = kepler_mean_to_ecc_anomaly (manom, e, &info);

        // Did solve fail?
        if (info != 0) {
            if (!status_) status_ = info;
            return;
        }

        cpsi = cos(psi);
        spsi = sin(psi);

        // Compute the true anomaly.
        d = 1.0 - e * cpsi;
        if (d == 0.0){
            if (!status_) status_ = 3;
            return;
        }
        cth = (cpsi - e) / d;
        sth = sqrt(1 - cth * cth);

        // Compute the radius and derivative.
        denom = 1.0 / (1 + e * cth);
        r = a * (1 - e * e) * denom;

        // Compute the coordinates in the plane of the orbit.
        sgn = (spsi >= 0) - (spsi < 0);
        x = r * cth;
        y = r * sth * sgn;

        // Rotate by pomega.
        xp =  x * cpomega + y * spomega;
        yp = -x * spomega + y * cpomega;

        // Rotate by the inclination angle.
        xsx = xp * six;

        // Compute the positions.
        pos[0] = xp * cix;
        pos[1] = yp * ciy - xsx * siy;
        pos[2] = yp * siy + xsx * ciy;
    };

    virtual double operator () (double t) {
        int i, l = mp_.size();
        double z, lam = 1.0, pos[3] = {0, 0, 0};
        for (i = 0; i < l; ++i) {
            // Solve Kepler's equation for the position of the body.
            solve(t, i, pos);

            // Find the impact parameter in the plane of the sky.
            z = sqrt(pos[1]*pos[1] + pos[2]*pos[2]);

            // Is the body transiting or occulting.
            if (pos[0] > 0.0) {
                // If transiting, use the stellar limb darkening profile.
                lam *= ld_(ror_[i], z/rstar_);
            } else if (occ_[i] > 0.0) {
                // If occulting, use the planet's limb darkening profile.
                double l = (*(pld_[i]))(iror_[i], z/r_[i]);
                lam *= 1.0 + occ_[i] * (l - 1);
            }

            // // Include Doppler beaming, ellipsoidal variations, and reflection.
            // if (ae_[i] > 0.0 || ab_[i] > 0.0 || ar_[i] > 0.0 || fabs(zp_[i]) > 0.0) {
            //     double phi = fmod(t-t0_[i], periods_[i])/periods_[i], f = zp_[i];
            //     if (ae_[i] > 0.0)
            //         f -= ae_[i] * cos(4*M_PI*phi - lag_[i]);
            //     if (ab_[i] > 0.0)
            //         f += ab_[i] * sin(2*M_PI*phi);
            //     if (ar_[i] > 0.0) {
            //         double cz = -cix_[i] * sin(2*M_PI*phi), z0 = acos(cz);
            //         f -= ar_[i] * (sin(z0) + (M_PI - z0) * cz) / M_PI;
            //     }
            //     lam *= 1+f;
            // }
        }
        return lam;
    };

protected:

    int status_;
    L ld_;
    double mstar_, rstar_;
    std::vector<LimbDarkening*> pld_;
    std::vector<double> occ_, mp_, r_, ror_, iror_, a_, t0_, e_, pomega_, ix_,
                        iy_, periods_, dmanomdt_, t1s_, cpom_, spom_,
                        cix_, six_, ciy_, siy_;
                        // zp_, ae_, lag_, ab_, ar_;

};
};

#endif
