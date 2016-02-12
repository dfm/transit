#ifndef _TRANSIT_KEPLER_H_
#define _TRANSIT_KEPLER_H_

#include <cfloat>
#include <vector>
#include <cstddef>

#define KEPLER_MAX_ITER 200
#define KEPLER_CONV_TOL 1.48e-7

using std::vector;

//
// Newton's constant in $R_\odot^3 M_\odot^{-1} {days}^{-2}$.
//
#define G_GRAV 2945.462538537765

namespace transit {

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
        if (fabs((psi - psi0) / psi) <= KEPLER_CONV_TOL) return psi;

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
int solve_kepler (const T& M,
                  const T& e, const T& cos_omega, const T& sin_omega,
                  const T& a, const T& cos_i, const T& sin_i,
                  T* pos)
{
    // Solve Kepler's equation.
    int info;
    T E = mean_to_ecc_anomaly (M, e, &info);
    if (info) return info;

    // Convert to Cartesian coordinates.
    T r = a * (1.0 - e * cos(E)),
      f = 2.0 * atan2(sqrt(1.0 + e) * tan(0.5 * E), sqrt(1.0 - e)),
      cf = cos(f), sf = sin(f),
      swpf = sin_omega * cf + cos_omega * sf,
      cwpf = cos_omega * cf - sin_omega * sf;
    pos[0] = r * cwpf;
    pos[1] = r * swpf * cos_i;
    pos[2] = r * swpf * sin_i;
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
    void position (const T* const params, const double t, T* pos) {
        // Access the parameters.
        T tref               = params[1],
          e                  = params[2],
          cos_omega          = params[3],
          sin_omega          = params[4],
          a                  = params[5],
          cos_i              = params[6],
          sin_i              = params[7],
          two_pi_over_period = params[8],
          M                  = (t - tref) * two_pi_over_period;

        status_ = solve_kepler (M, e, cos_omega, sin_omega, a, cos_i, sin_i, pos);
    };

    template <typename T>
    T* reparameterize (const T* const params) const {
        //
        // params:
        //
        //   [ln_fstar, ln_rstar, ln_mstar] +
        //   [ln_r, ln_m, t0, sqrt(e)*cos(pomega), sqrt(e)*sin(pomega),
        //    sqrt(a)*cos(ix), sqrt(a)*sin(ix)] * N_BODY +
        //   [q_1, q_2]
        //
        //  (TODO: [ln(q_1 / (1 - q_1)), ln(q_2 / (1 - q_2))])
        //
        T* result = new T[5 + 9 * n_body_];
        T m_central = exp(params[2]);

        // Central parameters.
        result[0]              = exp(params[0]);                             // Flux
        result[1]              = exp(params[1]);                             // Radius
        result[2+9*n_body_]    = 1.0 / (1.0 + exp(-params[3+7*n_body_]));    // q1
        result[2+9*n_body_+1]  = 1.0 / (1.0 + exp(-params[3+7*n_body_+1]));  // q2
        result[2+9*n_body_+2]  = 1.0 / (1.0 + exp(-params[3+7*n_body_+2]));  // dilution

        int i, j, n;
        for (i = 0, j = 3, n = 2; i < n_body_; ++i, j += 7, n += 9) {
            // Access the parameters.
            T mass        = exp(params[j + 1]),
              t0          = params[j + 2],
              sqrt_e_cosw = params[j + 3],
              sqrt_e_sinw = params[j + 4],
              sqrt_a_cosi = params[j + 5],
              sqrt_a_sini = params[j + 6],
              tref,
              e, sqrt_e,
              sin_omega, cos_omega,
              a, sqrt_a,
              two_pi_over_period,
              E0;

            // Compute the period.
            a = sqrt_a_cosi * sqrt_a_cosi + sqrt_a_sini * sqrt_a_sini;
            sqrt_a = sqrt(a);
            two_pi_over_period = sqrt((m_central + mass) * G_GRAV / (a*a*a));

            // Compute the eccentricity and special case circular orbits.
            e = sqrt_e_cosw * sqrt_e_cosw + sqrt_e_sinw * sqrt_e_sinw;
            if (e > DBL_EPSILON) {
                sqrt_e = sqrt(e);

                // omega
                sin_omega = sqrt_e_sinw / sqrt_e;
                cos_omega = sqrt_e_cosw / sqrt_e;

                // Compute the reference time using the time of transit.
                //  -> f = 0.5*pi - w
                //    -> sin(f) = cos(w)
                //    -> cos(f) = sin(w)
                //  -> tan(0.5*f) = sin(f) / (1 + cos(f))
                //                = cos(w) / (1 + sin(w))
                E0 = 2.0 * atan2(sqrt(1.0 - e) * cos_omega,
                                 sqrt(1.0 + e) * (1.0 + sin_omega));
                tref = t0 - (E0 - e * sin(E0)) / two_pi_over_period;
            } else {
                sin_omega = T(0.0);
                cos_omega = T(1.0);
                tref = t0 - T(0.5 * M_PI) / two_pi_over_period;
            }

            result[n]   = exp(params[j]);        // Radius
            result[n+1] = tref;                  // Reference time
            result[n+2] = e;                     // Eccentricity
            result[n+3] = cos_omega;             // ...
            result[n+4] = sin_omega;             // ...
            result[n+5] = a;                     // Semi-major axis
            result[n+6] = sqrt_a_cosi / sqrt_a;  // Inclination
            result[n+7] = sqrt_a_sini / sqrt_a;  // ...
            result[n+8] = two_pi_over_period;    // Time factor
        }

        return result;
    };

    template <typename T>
    T operator () (const T* const params, const double t) {
        int i, n, nld = 2 + 9 * n_body_;
        T z, lam = T(1.0), pos[3] = {T(0.0), T(0.0), T(0.0)};
        for (i = 0, n = 2; i < n_body_; ++i, n += 9) {
            // Solve Kepler's equation for the position of the body.
            position(&(params[n]), t, pos);

            // Find the impact parameter in the plane of the sky.
            z = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);

            // Is the body transiting?
            if (pos[2] > 0.0) {
                // If transiting, use the stellar limb darkening profile.
                lam *= ld_(&(params[nld]), params[n] / params[1], z / params[1]);
            }
        }
        return params[0] * ((1.0 - params[nld+2]) * lam + params[nld+2]);
    };

private:
    int status_, n_body_;
    L ld_;
};

};

#endif
