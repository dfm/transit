#ifndef _TRANSIT_KEPLER_H_
#define _TRANSIT_KEPLER_H_

#include <cfloat>
#include <vector>
#include <cstddef>

#define KEPLER_MAX_ITER 200
#define KEPLER_CONV_TOL 1.48e-10

using std::vector;

template <typename T> int sgn(T val) {
    return (T(0) <= val) - (val < T(0));
}

//
// Newton's constant in $R_\odot^3 M_\odot^{-1} {days}^{-2}$.
//
#define G_GRAV 2945.4625385377644

namespace transit {

template <typename T>
T ecc_to_mean_anomaly (const T& psi, const T& e) {
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

template <typename T, int N>
Jet<T, N> mean_to_ecc_anomaly (Jet<T, N> M, Jet<T, N> e, int *info) {
    T psi = mean_to_ecc_anomaly(M.a, e.a, info);
    return Jet<T, N>(
        psi,
        (M.v + e.v * sin(psi)) / (1.0 - e.a * cos(psi))
    );
}

template <typename T>
void ecc_anomaly_to_cartesian (const T& a, const T& e, const T& psi,
                               T* x, T* y)
{
    T hpsi = 0.5 * psi, shp = sin(hpsi), chp = cos(hpsi),
      theta = 2.0 * atan2(sqrt(1.0 - e) * chp, sqrt(1.0 + e) * shp),
      cth = cos(theta),
      factor = a * (e * e - 1.0) / (1.0 + e * cos(theta));
    *x = factor * cth;
    *y = factor * sin(theta);
}

template <typename T>
int solve_kepler (const T& manom,
                  const T& a, const T& e,
                  const T& sp, const T& cp,
                  const T& six, const T& cix,
                  const T& siy, const T& ciy,
                  T* pos)
{
    int info;
    T psi = mean_to_ecc_anomaly (manom, e, &info);

    // Did solve fail?
    if (info) return info;

    // Convert to Cartesian coordinates.
    T x, y;
    ecc_anomaly_to_cartesian (a, e, psi, &x, &y);

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

    KeplerSolver () : status_(0), n_body_(0), ld_() {};
    int get_status () const { return status_; };

    void set_n_body (const int n) { n_body_ = n; };
    int get_n_body () const { return n_body_; };

    template <typename T>
    void position (const T& m_central, const T* const params,
                   const double t, T* pos) {
        // Access the parameters.
        T a_body = params[2],
          e_body = params[4],
          cpom = params[5],
          spom = params[6],
          cix = params[7],
          six = params[8],
          ciy = params[9],
          siy = params[10],
          manom = t * params[11] + params[12];

        status_ = solve_kepler (manom, a_body, e_body, spom, cpom, six, cix, siy, ciy, pos);
    };

    template <typename T>
    T* reparameterize (const T* const params) const {
        //
        // params:
        //
        //   [ln_fstar, ln_rstar, ln_mstar] +
        //   [ln_r, ln_m, t0, sqrt(e)*cos(pomega), sqrt(e)*sin(pomega),
        //    sqrt(a)*cos(ix), sqrt(a)*sin(ix), iy] * N_BODY +
        //   [q_1, q_2]
        //
        //  (TODO: [ln(q_1 / (1 - q_1)), ln(q_2 / (1 - q_2))])
        //
        T* result = new T[5 + 13 * n_body_];
        T m_central = exp(params[2]);

        // Central parameters.
        result[0] = exp(params[0]);
        result[1] = exp(params[1]);
        result[2] = m_central;
        result[3+13*n_body_] = params[3+8*n_body_];
        result[3+13*n_body_+1] = params[3+8*n_body_+1];

        int i, j, n;
        for (i = 0, j = 3, n = 3; i < n_body_; ++i, j += 8, n += 13) {
            // Access the parameters.
            T m_body = exp(params[j + 1]),
              t0_body = params[j + 2],
              ecosp = params[j + 3],
              esinp = params[j + 4],
              acosi = params[j + 5],
              asini = params[j + 6],
              iy = params[j + 7],
              e_body, sqrte, a_body, sqrta,
              psi0, period_over_2pi;

            e_body = ecosp * ecosp + esinp * esinp;
            a_body = acosi*acosi + asini*asini;
            sqrta = sqrt(a_body);

            period_over_2pi = sqrt(a_body*a_body*a_body / (m_central+m_body) / G_GRAV);

            result[n]   = exp(params[j]);
            result[n+1] = m_body;
            result[n+2] = a_body;
            result[n+3] = t0_body;

            if (e_body > 0.0) {
                sqrte = sqrt(e_body);
                psi0 = 2.0 * atan2(-sqrt(1.0 - e_body) * esinp / (sqrte + ecosp),
                                   sqrt(1.0 + e_body));
                psi0 = psi0 - e_body * sin(psi0);
                result[n+4] = e_body;
                result[n+5]  = ecosp / sqrte;
                result[n+6]  = esinp / sqrte;
            } else {
                psi0 = T(0.0);
                result[n+4] = T(0.0);
                result[n+5]  = T(1.0);
                result[n+6]  = T(0.0);
            }

            result[n+7]  = acosi / sqrta;
            result[n+8]  = asini / sqrta;
            result[n+9]  = cos(iy);
            result[n+10] = sin(iy);

            result[n+11] = 1.0 / period_over_2pi;
            result[n+12] = -t0_body / period_over_2pi + psi0;
        }

        return result;
    };

    template <typename T>
    T operator () (const T* const params, const double t) {
        int i, n, nld = 3 + 13 * n_body_;
        T z, lam = params[0], pos[3] = {T(0.0), T(0.0), T(0.0)};
        for (i = 0, n = 3; i < n_body_; ++i, n += 13) {
            // Solve Kepler's equation for the position of the body.
            position(params[2], &(params[n]), t, pos);

            // Find the impact parameter in the plane of the sky.
            z = sqrt(pos[1]*pos[1] + pos[2]*pos[2]);

            // Is the body transiting?
            if (pos[0] > 0.0) {
                // If transiting, use the stellar limb darkening profile.
                lam *= ld_(&(params[nld]), params[n] / params[1], z / params[1]);
            }
        }
        return lam;
    };

private:
    int status_, n_body_;
    L ld_;
};

};

#endif
