#ifndef _KEPLER_H_
#define _KEPLER_H_

#include <cmath>
#include <cfloat>
#include <vector>

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

    int add_body (double m, double r, double a, double t0, double e,
                  double pomega, double ix, double iy)
    {
        if (e < 0.0 || fabs(e - 1.0) <= DBL_EPSILON) return 1;

        mp_.push_back(m);
        ror_.push_back(r/rstar_);
        a_.push_back(a);
        t0_.push_back(t0);
        e_.push_back(e);
        pomega_.push_back(pomega);
        ix_.push_back(ix);
        iy_.push_back(iy);

        // Pre-compute some constant factors.
        double period = 2 * M_PI * sqrt(a*a*a/G_GRAV/(mstar_+m)),
               psi0 = 2 * atan2(tan(0.5 * pomega), sqrt((1 + e) / (1 - e)));
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
               // drdt, dthdt, dxdt, dydt, dxpdt, dypdt,
               // dxsxdt;

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
        // denom = sqrt(1 - e * e);
        // drdt = dmanomdt * a * e * sth / denom;
        // dthdt = dmanomdt * a * a * denom / r / r;

        // Compute the coordinates in the plane of the orbit.
        sgn = (spsi >= 0) - (spsi < 0);
        x = r * cth;
        y = r * sth * sgn;
        // dxdt = drdt * cth - r * sth * dthdt;
        // dydt = (drdt * sth + r * cth * dthdt) * sgn;

        // Rotate by pomega.
        xp =  x * cpomega + y * spomega;
        yp = -x * spomega + y * cpomega;

        // Rotate by the inclination angle.
        xsx = xp * six;

        // Compute the positions.
        pos[0] = xp * cix;
        pos[1] = yp * ciy - xsx * siy;
        pos[2] = yp * siy + xsx * ciy;

        // // Rotate the velocities.
        // dxpdt =  dxdt * cpomega + dydt * spomega;
        // dypdt = -dxdt * spomega + dydt * cpomega;
        // dxsxdt = dxpdt * six;
        // // Compute the velocities.
        // coords[6*i+3] = sgn * (dxpdt * cix);
        // coords[6*i+4] = sgn * (dypdt * ciy - dxsxdt * siy);
        // coords[6*i+5] = sgn * (dypdt * siy + dxsxdt * ciy);
    };

    double operator () (double t) {
        int i, l = mp_.size();
        double lam = 1.0, pos[3];
        for (i = 0; i < l; ++i) {
            solve(t, i, pos);
            if (pos[0] > 0.0)
                lam *= ld_(ror_[i], sqrt(pos[1]*pos[1] + pos[2]*pos[2]));
        }
        return lam;
    };

private:

    L ld_;
    int status_;
    double mstar_, rstar_;
    std::vector<double> mp_, ror_, a_, t0_, e_, pomega_, ix_, iy_,
                        periods_, dmanomdt_, psi0s_, t1s_, cpom_, spom_,
                        cix_, six_, ciy_, siy_;

};
}

#endif
