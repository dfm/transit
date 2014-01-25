#ifndef _DRIVER_H_
#define _DRIVER_H_

int ldlc (double u1, double u2, double mstar, double rstar,
          int nplanets, double *m, double *r, double *a, double *t0,
          double *e, double *pomega, double *ix, double *iy,
          double texp, double tol, int maxdepth,
          int N, double *t, double *lam);

#endif
