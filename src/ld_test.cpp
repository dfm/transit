#include <iostream>

#include "limb_darkening.h"
#include "numerical_limb_darkening.h"

using transit::QuadraticLaw;
using transit::NumericalLimbDarkening;
using transit::QuadraticLimbDarkening;

int main ()
{
    double q1 = 0.9999, q2 = 0.9999,
           u1 = 2*q1*q2, u2 = q1*(1-2*q2);

    QuadraticLaw* law = new QuadraticLaw(u1, u2);
    NumericalLimbDarkening ld(law, 1e-6, 1000, 100);
    QuadraticLimbDarkening gld(u1, u2);

    int n = 0;
    double p = 0.1, z, norm = 0.0, mx = -INFINITY;

    for (z = 0.0; z < 1.1+p; z += 1.438956e-4, ++n) {
        double v, v0;
        v = ld(p, z);
        v0 = gld(p, z);
        if (v0 > 1) std::cout << z - (1+p) << " " << v0 - 1 << std::endl;
        norm += (v - v0) * (v - v0);
        if (fabs(v-v0) > mx) mx = fabs(v - v0);
    }

    std::cout << norm / n << std::endl;
    std::cout << mx << std::endl;

    delete law;
    return 0;
}
