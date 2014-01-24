#include <cstdio>
#include "driver.h"

#define N 50000


int main()
{
    int i;
    double m[] = {0.0, 0.0},
           r[] = {0.01, 0.03},
           a[] = {100.0, 50.0},
           t0[] = {1.0, 2.0},
           e[] = {0.0, 0.99},
           pomega[] = {0.0, 1.5},
           ix[] = {0.0, 0.0},
           iy[] = {0.0, 0.0},
           t[N], lam[N];

    for (i = 0; i < N; ++i) t[i] = i*0.001;
    int info = ldlc (0.3, 0.1, 1.0, 1.0, 2, m, r, a, t0, e, pomega, ix, iy, 1e-2, 0.1, 2, N, t, lam);
    /* printf("info = %d\n", info); */
    /* return 0; */
    for (i = 0; i < N; ++i) printf("%e %e\n", t[i], lam[i]);

    return 0;
}
