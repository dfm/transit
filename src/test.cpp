#include <cstdio>

#include "quad.h"
#include "solver.h"
#include "integrator.h"

using transit::Integrator;
using transit::SimpleSolver;
using transit::QuadraticLimbDarkening;

int main()
{
    double p = 0.01, t;
    SimpleSolver solver(10., 5., 0.1, p, 0.5);
    QuadraticLimbDarkening ld(0.3, 0.2);

    Integrator<SimpleSolver, QuadraticLimbDarkening> integrator(solver, ld, 1e-1, 10);
    for (t = 0.0; t < 10; t += 0.001) {
        printf("%20.16e %20.16e\n", t, integrator(t, 1e-2));
    }
    return 0;
}
