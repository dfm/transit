#ifndef _TRANSIT_SIMPLE_H_
#define _TRANSIT_SIMPLE_H_

#include <cmath>
#include "ceres/jet.h"

using ceres::Jet;

namespace transit {

inline double wrap (double x, double rng) { return x - rng * floor(x / rng); }
template <typename T, int N>
Jet<T, N> wrap (const Jet<T, N>& x, const Jet<T, N>& rng) {
    return Jet<T, N>(wrap(x.a, rng.a), x.v);
}

template <class L>
class SimpleSolver {

public:

    SimpleSolver () {};
    int get_status () const { return 0; };
    L& get_limb_darkening () { return ld_; };

    template <typename T>
    T operator () (const T* const params, const double t) const {
        T ror = params[0],
          one_plus_ror = 1.0 + ror,
          period = params[1],
          t0 = params[2],
          half_period = 0.5 * period,
          b2 = params[3] * params[3],
          duration = params[4],
          factor = 2.0 * sqrt(one_plus_ror * one_plus_ror - b2)/duration,
          x = factor * (wrap(t+half_period - t0, period) - half_period);
        return ld_(&(params[5]), ror, sqrt(b2 + x*x));
    };

private:
    L ld_;
};

};

#endif
