#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

namespace transit {
template <class S>
class Integrator {
public:
    Integrator (S solver, double tol, int maxdepth)
        : tol_(tol), maxdepth_(maxdepth), solver_(solver) {};

    int get_status () const { return solver_.get_status(); };

    double integrate (double f0, double t, double texp, int depth) {
        double st = texp/3.0,
               tp = t+st,
               tm = t-st,
               fp = solver_(tp),
               fm = solver_(tm),
               d = fabs((fp - 2*f0 + fm) / (fp - fm));

        if (d > tol_ && depth < maxdepth_) {
            fp = integrate(fp, tp, st, depth+1);
            f0 = integrate(f0, t, st, depth+1);
            fm = integrate(fm, tm, st, depth+1);
        }
        return (f0+fp+fm) / 3.0;
    };

    double operator () (double t, double texp) {
        return integrate(solver_(t), t, texp, 0);
    };

private:
    double tol_;
    int maxdepth_;
    S solver_;
};
}

#endif
