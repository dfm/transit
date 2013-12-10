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
    int i, j;
    double d, tmp, delta = 0.0, p, mu1 = 2*sqrt(q1)*q2,
           mu2 = sqrt(q1)*(1-2*q2),
           *z0 = (double*)malloc(n * sizeof(double)),
           *muo1 = (double*)malloc(n * sizeof(double)),
           *mu0 = (double*)malloc(n * sizeof(double));

    for (i = 0; i < n; ++i) z0[i] = 2.0 * (double(i)+M_PI/3)/n;
    for (i = 0; i < np; ++i) {
        p = sqrt((double(i)+1)/np);
        occultquad_ (z0, &mu1, &mu2, &p, muo1, mu0, &n);
        for (j = 0; j < n; ++j) {
            tmp = ldlc(p, z0[j], mu1, mu2);
            d = tmp - muo1[j];
            delta += d*d;
        }
    }

    free(z0);
    free(muo1);
    free(mu0);

    return delta;
}

int main()
{
    printf("%e\n", test_quad(0.5, 0.2, 500, 100));
    return 0;
}
