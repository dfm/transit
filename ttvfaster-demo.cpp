// compile using something like:
//
//  g++ ttvfaster-demo.cpp -o ttvfaster-demo -Itransit/include -I/usr/local/include/eigen3
//
// you'll need the Eigen headers somewhere.


#include <iostream>
#include "ceres/jet.h"
#include "ttvfaster.h"

using ceres::Jet;
using namespace transit::ttvfaster;

int main () {
    double params[15] = {
        1.0,
        0.00001027,  66.03300476, -0.00654273, 1.57079637, -1.57079637, -0.05891280, 142.63816833,
        0.00041269, 125.85229492, -0.00096537, 1.57079637, -1.57079637, -0.00953292, 236.66268921
    };

    // Start with the standard result.
    TTVFasterResult<double>* result = TTVFaster<double>().compute_times(2, params, 0.0, 1600.0, 6);
    for (unsigned i = 0; i < 2; ++i) {
        for (unsigned j = 0; j < result->get_ntransits(i); ++j) {
            std::cout << i << " " << j << " " << result->get_times(i)[j] << std::endl;
        }
    }
    delete result;

    // Then the version with the gradients.
    Jet<double, 15>* params2 = new Jet<double, 15>[15];
    for (unsigned i = 0; i < 15; ++i) {
        params2[i] = Jet<double, 15>(params[i], i);
    }
    TTVFasterResult<Jet<double, 15> >* result2 =
        TTVFaster<Jet<double, 15> >().compute_times(2, params2, 0.0, 1600.0, 6);
    for (unsigned i = 0; i < 2; ++i) {
        for (unsigned j = 0; j < result2->get_ntransits(i); ++j) {
            std::cout << i << " " << j << " " << result2->get_times(i)[j] << std::endl;
        }
    }
    delete result2;

    return 0;
}
