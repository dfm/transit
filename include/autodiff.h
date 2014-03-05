#ifndef _AUTODIFF_H_
#define _AUTODIFF_H_

#include <cmath>
#include <Eigen/Core>

namespace transit {

template <int N>
struct AutoDiff {

    AutoDiff () : a() {
        v.setZero();
    };
    explicit AutoDiff (const double& value) {
        a = value;
        v.setZero();
    };
    template<typename Derived>
    AutoDiff (const double& value, const Eigen::DenseBase<Derived>& vIn)
        : a(value), v(vIn) {}
    AutoDiff (const double& value, int k) {
        a = value;
        v.setZero();
        v[k] = 1.0;
    }

    // The scalar and vector parts.
    double a;
    Eigen::Matrix<double, N, 1, Eigen::DontAlign> v;

};

// Unary operators.
template<int N> inline
AutoDiff<N> const& operator+(const AutoDiff<N>& f) {
    return f;
}
template<int N> inline
AutoDiff<N> operator-(const AutoDiff<N>& f) {
    return AutoDiff<N>(-f.a, -f.v);
}

// Overload standard binary operations.

// Add.
template<int N> inline
AutoDiff<N> operator+(const AutoDiff<N>& f, const AutoDiff<N>& g) {
    return AutoDiff<N>(f.a + g.a, f.v + g.v);
}
template<int N> inline
AutoDiff<N> operator+(const AutoDiff<N>& f, double s) {
    return AutoDiff<N>(f.a + s, f.v);
}
template<int N> inline
AutoDiff<N> operator+(double s, const AutoDiff<N>& f) {
    return AutoDiff<N>(f.a + s, f.v);
}

// Subtract.
template<int N> inline
AutoDiff<N> operator-(const AutoDiff<N>& f, const AutoDiff<N>& g) {
    return AutoDiff<N>(f.a - g.a, f.v - g.v);
}
template<int N> inline
AutoDiff<N> operator-(const AutoDiff<N>& f, double s) {
    return AutoDiff<N>(f.a - s, f.v);
}
template<int N> inline
AutoDiff<N> operator-(double s, const AutoDiff<N>& f) {
    return AutoDiff<N>(s - f.a, f.v);
}

// Product.
template<int N> inline
AutoDiff<N> operator*(const AutoDiff<N>& f, const AutoDiff<N>& g) {
    return AutoDiff<N>(f.a * g.a, f.a * g.v + f.v * g.a);
}
template<int N> inline
AutoDiff<N> operator*(const AutoDiff<N>& f, double s) {
    return AutoDiff<N>(f.a * s, f.v * s);
}
template<int N> inline
AutoDiff<N> operator*(double s, const AutoDiff<N>& f) {
    return AutoDiff<N>(s * f.a, f.v * s);
}

// Division.
template<int N> inline
AutoDiff<N> operator/(const AutoDiff<N>& f, const AutoDiff<N>& g) {
    // This uses:
    //
    //   a + u   (a + u)(b - v)   (a + u)(b - v)
    //   ----- = -------------- = --------------
    //   b + v   (b + v)(b - v)        b^2
    //
    // which holds because v*v = 0.
    const double g_a_inverse = 1.0 / g.a;
    const double f_a_by_g_a = f.a * g_a_inverse;
    return AutoDiff<N>(f.a * g_a_inverse,
                       (f.v - f_a_by_g_a * g.v) * g_a_inverse);
}
template<int N> inline
AutoDiff<N> operator/(double s, const AutoDiff<N>& g) {
    const double minus_s_g_a_inverse2 = -s / (g.a * g.a);
    return AutoDiff<N>(s / g.a, g.v * minus_s_g_a_inverse2);
}
template<int N> inline
AutoDiff<N> operator/(const AutoDiff<N>& f, double s) {
    const double s_inverse = 1.0 / s;
    return AutoDiff<N>(f.a * s_inverse, f.v * s_inverse);
}

// Binary comparison operators.
#define TRANSIT_DEFINE_COMPARISON_OPERATOR(op) \
template<int N> inline \
bool operator op(const AutoDiff<N>& f, const AutoDiff<N>& g) { \
    return f.a op g.a; \
} \
template<int N> inline \
bool operator op(const double& s, const AutoDiff<N>& g) { \
    return s op g.a; \
} \
template<int N> inline \
bool operator op(const AutoDiff<N>& f, const double& s) { \
    return f.a op s; \
}
TRANSIT_DEFINE_COMPARISON_OPERATOR( <  )  // NOLINT
TRANSIT_DEFINE_COMPARISON_OPERATOR( <= )  // NOLINT
TRANSIT_DEFINE_COMPARISON_OPERATOR( >  )  // NOLINT
TRANSIT_DEFINE_COMPARISON_OPERATOR( >= )  // NOLINT
TRANSIT_DEFINE_COMPARISON_OPERATOR( == )  // NOLINT
TRANSIT_DEFINE_COMPARISON_OPERATOR( != )  // NOLINT
#undef TRANSIT_DEFINE_COMPARISON_OPERATOR

// sqrt(a + h) ~= sqrt(a) + h / (2 sqrt(a))
inline double sqrt (double x) { return std::sqrt(x); }
template <int N> inline
AutoDiff<N> sqrt(const AutoDiff<N>& f) {
    const double tmp = sqrt(f.a);
    const double two_a_inverse = 1.0 / (2.0 * tmp);
    return AutoDiff<N>(tmp, f.v * two_a_inverse);
}

// acos(a + h) ~= acos(a) - 1 / sqrt(1 - a^2) h
inline double acos (double x) { return std::acos(x); }
template <int N> inline
AutoDiff<N> acos(const AutoDiff<N>& f) {
    const double tmp = -1.0 / sqrt(1.0 - f.a * f.a);
    return AutoDiff<N>(acos(f.a), tmp * f.v);
}

template <int N>
inline std::ostream &operator<<(std::ostream &s, const AutoDiff<N>& z) {
    return s << "[" << z.a << " ; " << z.v.transpose() << "]";
}

};

#endif
