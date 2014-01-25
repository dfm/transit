#ifndef _DRIVER_H_
#define _DRIVER_H_

#ifdef __cplusplus
extern "C" {
#endif

int ldlc (double u1, double u2, double mstar, double rstar,
          int nplanets, double *m, double *r, double *a, double *t0,
          double *e, double *pomega, double *ix, double *iy,
          double texp, double tol, int maxdepth,
          int N, double *t, double *lam);

#ifdef __cplusplus
}
#endif

#endif
