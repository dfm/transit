#include "solver.h"
#include "kepler.h"
#include "driver.h"
#include "integrator.h"
#include "limb_darkening.h"

using transit::Integrator;
using transit::KeplerSolver;
using transit::SimpleSolver;

using transit::LimbDarkening;
using transit::GeometricLimbDarkening;
using transit::QuadraticLimbDarkening;

int ldlc_kepler (double u1, double u2, double mstar, double rstar,
                 int nplanets, double *occ, double *m, double *r, double *a,
                 double *t0, double *e, double *pomega, double *ix, double *iy,
                 double texp, double tol, int maxdepth,
                 int N, double *t, double *lam)
{
    int i, info;

    GeometricLimbDarkening pld;
    QuadraticLimbDarkening ld (u1, u2);
    KeplerSolver<QuadraticLimbDarkening> solver(ld, mstar, rstar);

    // Add the planets to the system.
    for (i = 0; i < nplanets; ++i) {
        info = solver.add_body(&pld, occ[i], m[i], r[i], a[i], t0[i], e[i],
                               pomega[i], ix[i], iy[i]);
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

int ldlc_simple (double u1, double u2, double p, double t0, double tau,
                 double ror, double b,
                 double texp, double tol, int maxdepth,
                 int N, double *t, double *lam)
{
    QuadraticLimbDarkening ld (u1, u2);
    SimpleSolver<QuadraticLimbDarkening> solver(ld, p, t0, tau, ror, b);
    Integrator<SimpleSolver<QuadraticLimbDarkening> > integrator(solver, tol, maxdepth);

    // Compute the integrated light curve.
    int i;
    for (i = 0; i < N; ++i) {
        lam[i] = integrator(t[i], texp);
        if (integrator.get_status()) return integrator.get_status();
    }

    return 0;
}
