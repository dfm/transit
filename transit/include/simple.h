#ifndef _TRANSIT_SIMPLE_H_
#define _TRANSIT_SIMPLE_H_

#include <cmath>
#include "ceres/jet.h"

using ceres::Jet;

namespace transit {

inline double wrap (double x, double rng) { return x - rng * floor(x / rng); }

template <typename T>
T wrap (const T x, const T rng) {
    T y = x;
    while (y >= rng) y = y - rng;
    while (y < T(0.0)) y = y + rng;
    return y;
}

template <class L>
class SimpleSolver {

public:

    SimpleSolver () {};
    int get_status () const { return 0; };
    L& get_limb_darkening () { return ld_; };

    template <typename T>
    T* reparameterize (const T* const params) const {
        T* result = new T[7];
        T ror = exp(params[0]), rorp1 = 1.0 + ror;
        result[0] = ror;
        result[1] = exp(params[1]);
        result[2] = params[2];
        result[3] = params[3] * params[3];
        result[4] = 1.0 / (1.0 + exp(-params[5]));
        result[5] = 1.0 / (1.0 + exp(-params[6]));
        result[6] = 2.0 * sqrt(rorp1 * rorp1 - result[3]) / exp(params[4]);
        return result;
    }

    template <typename T>
    T operator () (const T* const params, const double t) const {
        T ror = params[0],
          period = params[1],
          t0 = params[2],
          b2 = params[3],
          factor = params[6],
          half_period = 0.5 * period,
          x = factor * (wrap(t + half_period - t0, period) - half_period);
        return ld_(&(params[4]), ror, sqrt(b2 + x*x));
    };

private:
    L ld_;
};

};

#endif
