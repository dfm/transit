#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <cmath>

namespace transit {
class SimpleSolver {

public:

    SimpleSolver (double period, double t0, double duration, double ror,
                  double impact)
        : period_(period), t0_(t0), duration_(duration),
          ror_(ror), impact_(impact)
    {
        opr_ = 1.0 + ror_;
        hp_ = 0.5 * period_;
        b2_ = impact_ * impact_;
        factor_ = 2.0 * sqrt(opr_*opr_ - b2_)/duration_;
    };

    double get_ror () const { return ror_; };

    double operator () (double t) const
    {
        double x = factor_ * (fmod(t+hp_-t0_, period_)-hp_);
        return sqrt(b2_ + x*x);
    };

private:

    double period_, t0_, duration_, ror_, impact_;
    double opr_, b2_, hp_, factor_;

};
}

#endif
