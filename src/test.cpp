#include "quad.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>

#ifdef __cplusplus
extern "C" {
#endif
extern void occultquad_ (double *z0, double *u1, double *u2, double *p,
                         double *muo1, double *mu0, int *nz);
#ifdef __cplusplus
}
#endif

double test_quad (double q1, double q2, int np, int n)
{
    int i, j, one = 1;
    double d, tmp, delta = 0.0, p, mu1 = 2*sqrt(q1)*q2,
           mu2 = sqrt(q1)*(1-2*q2), z0, muo1 = 0.0, mu0;

    for (i = 0; i < np; ++i) {
        p = sqrt((double(i)+1)/np);
        for (j = 0; j < n; ++j) {
            z0 = 2.0 * (double(j)+M_PI/3)/n;
            if (z0 >= p || z0 != 1-p) {
                tmp = ldlc(p, z0, mu1, mu2);
                occultquad_ (&z0, &mu1, &mu2, &p, &muo1, &mu0, &one);
                d = tmp - muo1;
                if (fabs(d) > 1e-5)
                    printf("%e %e %e\n", tmp, muo1, d);
                delta += d*d;
            }
        }
    }

    return sqrt(delta) / np / n;
}

int main()
{
    printf("%e\n", test_quad(0.5, 0.2, 100, 100));
    return 0;
}
