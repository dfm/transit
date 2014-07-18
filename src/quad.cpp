//
// Code adapted from A. Pal's publicly available implementation of the
// results from: Pal, A. 2008, MNRAS, 390, 281
//

#include <cmath>
#include <cfloat>
#include <boost/math/special_functions/ellint_3.hpp>

#include "limb_darkening.h"

using boost::math::ellint_1;
using boost::math::ellint_2;
using boost::math::ellint_3;


// Constants from Table A1.
#define PAL_AB  double a = (p-z)*(p-z), \
                       b = (p+z)*(p+z);
#define PAL_K   double k0=acos((p2+z2-1)/(2*p*z)), \
                       k1=acos((1+z2-p2)/(2*z));

#define PAL_CI  double ci = 2.0/(9.0*M_PI*sqrt(1.0-a)), \
                       cik = (1.0-5*z2+p2+a*b), \
                       cie = (z2+7*p2-4)*(1.0-a), \
                       cip = -3.0*(p+z)/(p-z); \
                fk = ci*cik; fe = ci*cie; fp = ci*cip; \
                kk = ee = pp = 1; \
                k = sqrt(4.0*z*p/(1.0-a)); \
                n = (a-b)/a; \
                f0 = p2; \
                f2 = 0.5*f0*(f0+2*z2);

#define PAL_CG  PAL_K \
                double cg = 1.0/(9.0*M_PI*sqrt(p*z)), \
                       cgk = ((1-b)*(2*b+a-3)-3*(p+z)*(p-z)*(b-2)), \
                       cge = 4*p*z*(z2+7*p2-4), \
                       cgp = -3*(p+z)/(p-z); \
                fk = cg*cgk; fe = cg*cge; fp = cg*cgp; \
                kk = ee = pp = 1; \
                k = sqrt((1-a)/(4*p*z)); \
                n = (a-1)/a; \
                f0 = (p2*k0+k1-sqrt(z2-0.25*(1+z2-p2)*(1+z2-p2)))/M_PI; \
                f2 = (k1+p2*(p2+2*z2)*k0-0.25*(1+5*p2+z2)*sqrt((1-a)*(b-1)))/(2.0*M_PI);


double transit::QuadraticLimbDarkening::operator () (double p, double z) const
{
    int kk = 0, ee = 0, pp = 0;
    double w0, w1, w2,
           f0 = 0, f1 = 0, fk = 0, fe = 0, fp = 0, f2 = 0, n = 0, k = 0,
           df, z2, p2;

    // Pre-compute the constant coefficients.
    w0 = (6-6*u1_-12*u2_) / (6-2*u1_-u2_);
    w1 = ( 6*u1_+12*u2_ ) / (6-2*u1_-u2_);
    w2 = (   6*u2_     ) / (6-2*u1_-u2_);
    z = fabs(z);
    z2 = z * z;
    p2 = p * p;

    // Run through the cases and compute the coefficients from Table A1 & A2.
    if (z <= 0.0 && p <= 1.0) {
        // M&A 10, Pal A
        f0 = p2;
        f1 = 2.0/3.0*(1.0-(1.0-f0)*sqrt(1.0-f0));
        f2 = 0.5*f0*f0;
    } else if (z <= p-1.0) {
        // M&A 11, Pal A_G
        f0 = 1.0;
        f1 = 2.0/3.0;
        f2 = 0.5;
    } else if (z < p && z < 1.0 - p - DBL_EPSILON) {
        // M&A 9, Pal B
        PAL_AB PAL_CI
        f1 = 2.0/3.0;
    } else if (z < p && fabs(z-1.0+p) <= DBL_EPSILON) {
        // M&A -, Pal B_T
        f0 = p2;
        f1 = 2.0/(3.0*M_PI)*acos(1-2*p)-4.0/(9.0*M_PI)*(3+2*p-8*f0)*sqrt(p*(1-p));
        f2 = 0.5*f0*(f0+2*z2);
    } else if (z < p) {
        // M&A 8, Pal B_G
        PAL_AB PAL_CG
        f1 = 2.0/3.0;
    } else if (fabs(z-p) <= DBL_EPSILON && z < 1.0-p-DBL_EPSILON) {
        // M&A 5, Pal C
        double t = 2.0/(9.0*M_PI);
        f0 = p2;
        f1 = 1.0/3.0;
        fk = t*(1.0-4.0*f0); kk = 1;
        fe = t*4*(2.0*f0-1.0); ee = 1;
        f2 = 1.5*f0*f0;
        k = 2*p;
    } else if (fabs(z-p) <= DBL_EPSILON && fabs(z-1.0+p) <= DBL_EPSILON) {
        // M&A 6, Pal C_T
        f0 = 0.25;
        f1 = 1.0/3.0-4.0/(9.0*M_PI);
        f2 = 3.0/32.0;
    } else if (fabs(z-p) <= DBL_EPSILON) {
        // M&A 7, Pal C_G
        PAL_AB PAL_K
        f0 = (p2*k0+k1-sqrt(z2-0.25*(1+z2-p2)*(1+z2-p2)))/M_PI;
        f1 = 1.0/3.0;
        fk = -(1-4*p2)*(3-8*p2)/(9*M_PI*p); kk = 1;
        fe = 16*p*(2*p2-1)/(9*M_PI); ee = 1;
        f2 = (k1+p2*(p2+2*z2)*k0-0.25*(1+5*p2+z2)*sqrt((1-a)*(b-1)))/(2.0*M_PI);
        k = 1.0/(2.0*p);
    } else if (z < 1-p-DBL_EPSILON) {
        // M&A 3, Pal D
        PAL_AB PAL_CI
    } else if (fabs(z-1.0+p) <= DBL_EPSILON) {
        // M&A 4, Pal E
        f0 = p2;
        f1 = 2.0/(3.0*M_PI)*acos(1-2*p)-4.0/(9.0*M_PI)*(3+2*p-8*f0)*sqrt(p*(1-p));
        f2 = 0.5*f0*(f0+2*z2);
    } else if (z < 1+p-DBL_EPSILON) {
        // M&A 2, Pal F
        PAL_AB PAL_CG
    }

    // Compute the base flux change.
    df = w0*f0+w1*f1+w2*f2;

    // Add in the elliptic integral terms.
    if (kk && fk != 0.0)
        df += w1 * fk * ellint_1<double>(k);
    if (ee && fe != 0.0)
        df += w1 * fe * ellint_2<double>(k);
    if (pp && fp != 0.0)
        df += w1 * fp * ellint_3<double, double>(-k, n);

    return 1.0 - df;
}
