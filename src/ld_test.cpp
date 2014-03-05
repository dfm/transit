#include <ctime>
#include <iostream>

#include "limb_darkening.h"

using transit::AutoDiff;
/* using transit::QuadraticLaw; */
using transit::LimbDarkeningLaw;
using transit::NumericalLimbDarkening2;
using transit::GeometricLimbDarkening;
/* using transit::QuadraticLimbDarkening; */

int main ()
{
    clock_t start, end;
    LimbDarkeningLaw* law = new LimbDarkeningLaw();
    /* QuadraticLaw* law = new QuadraticLaw(0.4, 0.2); */
    NumericalLimbDarkening2 ld(law, 1e-8, 500, 5);
    GeometricLimbDarkening gld;
    /* QuadraticLimbDarkening gld(0.4, 0.2); */

    double p = 0.1, z;

    for (z = 0.54586; z < 1.1+p; z += 0.001) {
        double v, v0;
        start = clock();
        v = ld(p, z);
        std::cout << double(clock()-start)/double(CLOCKS_PER_SEC) << " ";
        start = clock();
        v0 = gld(p, z);
        std::cout << double(clock()-start)/double(CLOCKS_PER_SEC) << " ";
        std::cout << z << "\t" << v << "\t" << v0 << "\t" << v - v0  << std::endl;
        break;
    }

    delete law;
    return 0;
}
