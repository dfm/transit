#include "quad.h"
#include "kepler.h"
#include "integrator.h"

using transit::Integrator;
using transit::KeplerSolver;
using transit::QuadraticLimbDarkening;

int ldlc (double u1, double u2, double mstar, double rstar,
          int nplanets, double *m, double *r, double *a, double *t0,
          double *e, double *pomega, double *ix, double *iy,
          double texp, double tol, int maxdepth,
          int N, double *t, double *lam)
{
    int i, info;

    QuadraticLimbDarkening ld (u1, u2);
    KeplerSolver<QuadraticLimbDarkening> solver(ld, mstar, rstar);

    // Add the planets to the system.
    for (i = 0; i < nplanets; ++i) {
        info = solver.add_body(m[i], r[i], a[i], t0[i], e[i], pomega[i],
                               ix[i], iy[i]);
        if (info) return info;
    }
    Integrator<KeplerSolver<QuadraticLimbDarkening> > integrator(solver, tol, maxdepth);

    // Compute the integrated light curve.
    for (i = 0; i < N; ++i) {
        lam[i] = integrator(t[i], texp);
        if (integrator.get_status()) return integrator.get_status();
    }

    return 0;
}
