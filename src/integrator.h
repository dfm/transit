#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

namespace transit {
template <class S, class L>
class Integrator {
public:
    Integrator (S solver, L ld, double tol, int maxdepth)
        : solver_(solver), ld_(ld), tol_(tol), maxdepth_(maxdepth) {};

    double evaluate (double t) const
    {
        return ld_(solver_.get_ror(), solver_(t));
    };

    double integrate (double f0, double t, double texp, int depth) const
    {
        double st = texp/3.0,
               tp = t+st,
               tm = t-st,
               fp = evaluate(tp),
               fm = evaluate(tm),
               d = fabs((fp - 2*f0 + fm) / (fp - fm));

        if (d > tol_ && depth < maxdepth_) {
            fp = integrate(fp, tp, st, depth+1);
            f0 = integrate(f0, t, st, depth+1);
            fm = integrate(fm, tm, st, depth+1);
        }
        return (f0+fp+fm) / 3.0;
    };

    double operator () (double t, double texp) const
    {
        return integrate(evaluate(t), t, texp, 0);
    };

private:
    double tol_;
    int maxdepth_;
    S solver_;
    L ld_;
};
}

#endif
