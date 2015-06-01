#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include <cmath>

namespace transit {
template <class S>
class Integrator {
public:
    Integrator () {};
    Integrator (S* solver, double tol, int maxdepth)
        : tol_(tol), maxdepth_(maxdepth), solver_(solver) {};

    int get_status () const { return solver_->get_status(); };

    double integrate (double f1, double f2, double t1, double t2, int depth) {
        double tmid = 0.5 * (t1 + t2),
               fmid = (*solver_)(tmid),
               fapprox = 0.5*(f1+f2),
               d = fabs(fmid - fapprox);
        if (d > tol_ && depth < maxdepth_) {
            double a = integrate(f1, fmid, t1, tmid, depth+1),
                   b = integrate(fmid, f2, tmid, t2, depth+1);
            return a + b;
        }
        return fapprox * (t2 - t1);
    };

    double integrate (double t, double texp) {
        texp = fabs(texp);
        if (texp > 0.0) {
            double dt = 0.5*texp,
                   t1 = t - dt,
                   t2 = t + dt;
            return integrate((*solver_)(t1), (*solver_)(t2), t1, t2, 0)/(t2-t1);
        }
        return (*solver_)(t);
    };

private:
    double tol_;
    int maxdepth_;
    S* solver_;
};
}

#endif
