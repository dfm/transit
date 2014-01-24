#include <cstdio>

#include "quad.h"
#include "kepler.h"
#include "integrator.h"

using transit::Integrator;
using transit::KeplerSolver;
using transit::QuadraticLimbDarkening;

int main()
{
    double p = 0.01, t;
    QuadraticLimbDarkening ld(0.3, 0.2);
    KeplerSolver<QuadraticLimbDarkening> solver(ld, 1.0, 1.0);
    solver.add_body(0.0, 0.01, 100., 5.0, 0.0, 0.0, 0.0, 0.0);
    solver.add_body(0.0, 0.03, 100., 6.0, 0.99, 0.0, 0.0, 0.0);

    Integrator<KeplerSolver<QuadraticLimbDarkening> > integrator(solver, 1e-1, 3);
    for (t = 0.0; t < 10; t += 0.001) {
        printf("%20.16e %20.16e\n", t, integrator(t, 1e-2));
    }
    return 0;
}
