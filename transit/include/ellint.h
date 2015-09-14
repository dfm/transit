#ifndef _TRANSIT_ELLINT_H_
#define _TRANSIT_ELLINT_H_

//
// Elliptic integrals computed following:
//  Bulirsch 1965, Numerische Mathematik, 7, 78
//  Bulirsch 1965, Numerische Mathematik, 7, 353
//
// And the implementation by E. Agol (private communication).
//

#include <cmath>
#include "ceres/jet.h"

using std::abs;
using ceres::Jet;

namespace transit {

#define ELLINT_CONV_TOL 1.0e-8
#define ELLINT_MAX_ITER 200

// K: 1.0 - k^2 >= 0.0
double ellint_1 (const double k) {
    unsigned i;
    double kc = sqrt(1.0 - k * k), m = 1.0, h;
    for (i = 0; i < ELLINT_MAX_ITER; ++i) {
        h = m;
        m += kc;
        if (abs(h - kc) / h <= ELLINT_CONV_TOL) break;
        kc = sqrt(h * kc);
        m *= 0.5;
    }
    if (i >= ELLINT_MAX_ITER)
        std::cerr << "Warning: ellint_1 did not converge; use result with caution\n";
    return M_PI / m;
}

// E: 1.0 - k^2 >= 0.0
double ellint_2 (const double k) {
    unsigned i;
    double b = 1.0 - k * k, kc = sqrt(b), m = 1.0, c = 1.0, a = b + 1.0, m0;
    for (i = 0; i < ELLINT_MAX_ITER; ++i) {
        b = 2.0 * (c * kc + b);
        c = a;
        m0 = m;
        m += kc;
        a += b / m;
        if (abs(m0 - kc) / m0 <= ELLINT_CONV_TOL) break;
        kc = 2.0 * sqrt(kc * m0);
    }
    if (i >= ELLINT_MAX_ITER)
        std::cerr << "Warning: ellint_2 did not converge; use result with caution\n";
    return M_PI_4 * a / m;
}

// Pi: 1.0 - k^2 >= 0.0 & 0.0 <= n < 1.0 (doesn't seem consistent for n < 0.0)
double ellint_3 (const double k, const double n) {
    unsigned i;
    double kc = sqrt(1.0 - k * k), p = sqrt(1.0 - n), m0 = 1.0, c = 1.0,
           d = 1.0 / p, e = kc, f, g;
    for (i = 0; i < ELLINT_MAX_ITER; ++i) {
        f = c;
        c += d / p;
        g = e / p;
        d = 2.0 * (f * g + d);
        p = g + p;
        g = m0;
        m0 = kc + m0;
        if (abs(1.0 - kc / g) <= ELLINT_CONV_TOL) break;
        kc = 2.0 * sqrt(e);
        e = kc * m0;
    }
    if (i >= ELLINT_MAX_ITER)
        std::cerr << "Warning: ellint_3 did not converge; use result with caution\n";
    return M_PI_2 * (c * m0 + d) / (m0 * (m0 + p));
}

// Gradients.
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

}; // namespace transit;

#endif  // _TRANSIT_ELLINT_H_
