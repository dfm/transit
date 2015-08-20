#ifndef _JET_EXT_H_
#define _JET_EXT_H_

#include <cmath>
#include "ceres/jet.h"
#include <boost/math/special_functions/ellint_3.hpp>

using ceres::Jet;

namespace transit {

inline double abs (double x) { return std::abs(x); }

template <typename T>
inline T _acos (const T& x) {
    if (x > T(1.0)) return acos(T(1.0));
    if (x < T(-1.0)) return acos(T(-1.0));
    return acos(x);
}

double ellint_1 (double z) { return boost::math::ellint_1<double>(z); }
double ellint_2 (double z) { return boost::math::ellint_2<double>(z); }
double ellint_3 (double z, double n) { return boost::math::ellint_3<double, double>(z, n); }

template <typename T, int N>
Jet<T, N> ellint_1 (const Jet<T, N>& z)
{
    T Kz = ellint_1(z.a),
      Ez = ellint_2(z.a),
      z2 = z.a*z.a;
    return Jet<T, N>(
        Kz,
        z.v * (Ez / (1.0 - z2) - Kz) / z.a
    );
}

template <typename T, int N>
Jet<T, N> ellint_2 (const Jet<T, N>& z) {
    T Kz = ellint_1(z.a),
      Ez = ellint_2(z.a);
    return Jet<T, N>(
        Ez,
        z.v * (Ez - Kz) / z.a
    );
}

template <typename T, int N>
Jet<T, N> ellint_3 (const Jet<T, N>& k, const Jet<T, N>& n) {
    T Kk = ellint_1(k.a),
      Ek = ellint_2(k.a),
      Pnk = ellint_3(k.a, n.a),
      k2 = k.a * k.a,
      n2 = n.a * n.a;
    return Jet<T, N>(
        Pnk,
        (n.v * 0.5*(Ek + (Kk*(k2-n.a) + Pnk*(n2-k2))/n.a) / (n.a-1.0) -
         k.v * k.a * (Ek / (k2 - 1.0) + Pnk)) / (k2-n.a)
    );
}

};

#endif
