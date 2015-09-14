#ifndef _TRANSIT_JET_EXT_H_
#define _TRANSIT_JET_EXT_H_

#include <cmath>

namespace transit {

inline double abs (double x) { return std::abs(x); }

template <typename T>
inline T _acos (const T& x) {
    if (x > T(1.0)) return acos(T(1.0));
    if (x < T(-1.0)) return acos(T(-1.0));
    return acos(x);
}

};

#endif
