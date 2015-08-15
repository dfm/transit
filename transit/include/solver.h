#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <cmath>

namespace transit {

inline double wrap (double x, double rng) { return x - rng * floor(x / rng); }

template <class L>
class SimpleSolver {

public:

    SimpleSolver () {};
    SimpleSolver (L* ld, double period, double t0, double duration, double ror,
                  double impact)
        : ld_(ld), period_(period), t0_(t0), duration_(duration),
          ror_(ror), impact_(impact)
    {
        opr_ = 1.0 + ror_;
        hp_ = 0.5 * period_;
        b2_ = impact_ * impact_;
        factor_ = 2.0 * sqrt(opr_*opr_ - b2_)/duration_;
    };
    ~SimpleSolver () {
        delete ld_;
    };

    int get_status () const { return 0; };

    double operator () (double t) const {
        double x = factor_ * (wrap(t+hp_-t0_, period_)-hp_);
        return (*ld_)(ror_, sqrt(b2_ + x*x));
    };

private:

    L* ld_;
    double period_, t0_, duration_, ror_, impact_;
    double opr_, b2_, hp_, factor_;

};
}

#endif
