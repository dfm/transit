#ifndef _TRANSIT_INTEGRATOR_H_
#define _TRANSIT_INTEGRATOR_H_

#include <cmath>
#include <vector>
#include "ceres/jet.h"

using std::vector;
using ceres::Jet;

namespace transit {

template <class S>
class Integrator {
public:
    Integrator ()
        : tol_(1.0e-2), maxdepth_(2), solver_() {};
    Integrator (double tol, int maxdepth)
        : tol_(tol), maxdepth_(maxdepth), solver_() {};

    int get_status () const { return solver_.get_status(); };
    S* get_solver () { return &solver_; };

    template <typename T>
    T integrate (const T* const params,
                 const T& f1, const T& f2,
                 double t1, double t2,
                 int depth) {
        double tmid = 0.5 * (t1 + t2);
        T fmid = solver_(params, tmid),
          fapprox = 0.5*(f1+f2),
          d = abs(fmid - fapprox);
        if (d > tol_ && depth < maxdepth_) {
            T a = integrate(params, f1, fmid, t1, tmid, depth+1),
              b = integrate(params, fmid, f2, tmid, t2, depth+1);
            return a + b;
        }
        return fapprox * (t2 - t1);
    };

    template <typename T>
    T integrate (const T* const params, double t, double texp) {
        texp = fabs(texp);
        if (texp > 0.0) {
            double dt = 0.5*texp,
                   t1 = t - dt,
                   t2 = t + dt;
            return integrate(params, solver_(params, t1), solver_(params, t2),
                             t1, t2, 0) / (t2 - t1);
        }
        return solver_(params, t);
    };

    int integrate (const double* params, int nt, const double* t,
                   double texp, double* lc) {
        int i, flag;
        double* reparams = solver_.reparameterize(params);

        for (i = 0; i < nt; ++i) {
            lc[i] = integrate(reparams, t[i], texp);
            flag = solver_.get_status();
            if (flag) return flag;
        }

        delete reparams;
        return 0;
    };

    template <int N>
    int gradient (const double* params, int nt, const double* t,
                  double texp, double* lc, double* gradient) {
        int i, j, flag = 0;
        Jet<double, N>* params2 = new Jet<double, N>[N];
        for (i = 0; i < N; ++i) params2[i] = Jet<double, N>(params[i], i);

        Jet<double, N>* reparams = solver_.reparameterize(params2);
        delete[] params2;

        for (i = 0; i < nt; ++i) {
            Jet<double, N> v = integrate(reparams, t[i], texp);
            flag = solver_.get_status();
            lc[i] = v.a;
            for (j = 0; j < N; ++j) gradient[i*N+j] = v.v[j];
            if (flag) break;
        }

        delete reparams;
        return flag;
    };

private:
    double tol_;
    int maxdepth_;
    S solver_;
};
};

#endif
