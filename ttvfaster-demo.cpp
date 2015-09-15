// compile using something like:
//
//  g++ ttvfaster-demo.cpp -o ttvfaster-demo -Itransit/include -I/usr/local/include/eigen3
//
// you'll need the Eigen headers somewhere.


#include <iostream>
#include <iomanip>
#include <cmath>
#include "ceres/jet.h"
#include "ttvfaster.h"

using ceres::Jet;
using namespace transit::ttvfaster;

int main () {
    double e1 = sqrt(pow(-0.00654273, 2) + pow(-0.05891280, 2)),
           omega1 = atan2(-0.05891280, -0.00654273),
           e2 = sqrt(pow(-0.00096537, 2) + pow(-0.00953292, 2)),
           omega2 = atan2(-0.00096537, -0.00953292);
    double params[15] = {
        log(1.0),
        log(0.00001027),  log(66.03300476), sqrt(e1)*cos(omega1), 1.57079637, -1.57079637, sqrt(e1)*sin(omega1), 142.63816833,
        log(0.00041269), log(125.85229492), sqrt(e2)*cos(omega2), 1.57079637, -1.57079637, sqrt(e2)*sin(omega2), 236.66268921
    };
    for (unsigned i = 0; i < 15; ++i)
        std::cout << std::setprecision(15) << params[i] << ", ";
    std::cout << std::endl;

    // Start with the standard result.
    TTVFaster<double> solver1;

    // Compute the number of transit.
    unsigned* n_transits = new unsigned[2],
            * starts = new unsigned[2];
    unsigned ntot = solver1.compute_ntransits(2, params, 0.0, 1600.0, n_transits, starts);

    // Run the solve.
    double* times = new double[ntot];
    solver1.compute_times(2, params, 0.0, 6, n_transits, starts, times);
    for (unsigned i = 0, n = 0; i < 2; ++i) {
        for (unsigned j = 0; j < n_transits[i]; ++j, ++n) {
            std::cout << i << " " << j << " " << times[n] << std::endl;
        }
    }

    // Now get the gradients.
    TTVFaster<Jet<double, 15> > solver2;
    Jet<double, 15>* params2 = new Jet<double, 15>[15];
    for (unsigned i = 0; i < 15; ++i) {
        params2[i] = Jet<double, 15>(params[i], i);
    }

    Jet<double, 15>* grad_times = new Jet<double, 15>[ntot];
    solver2.compute_times(2, params2, 0.0, 6, n_transits, starts, grad_times);
    for (unsigned i = 0, n = 0; i < 2; ++i) {
        for (unsigned j = 0; j < n_transits[i]; ++j, ++n) {
            std::cout << i << " " << j << " " << grad_times[n] << std::endl;
        }
    }

    delete grad_times;
    delete times;
    delete n_transits;
    delete starts;
    return 0;
}
