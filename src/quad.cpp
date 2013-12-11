#include <cmath>
#include <cfloat>
#include <boost/math/special_functions/ellint_3.hpp>
#include "quad.h"

/* #define VERBOSE */

using boost::math::ellint_1;
using boost::math::ellint_2;
using boost::math::ellint_3;

#define LAMBDA_1 { \
    double k = sqrt(0.25*(1-a)/(z*p)); \
    ld = ( ((1-b)*(2*b+a-3)-3*q*(b-2))*ellint_1<double>(k) \
           + 4*p*z*(z2+7*p2-4)*ellint_2<double>(k) \
           - 3*q/a*ellint_3<double, double>(-k, (a-1)/a) ) \
         / (9*M_PI*sqrt(p*z)); \
}

#define LAMBDA_2 { \
    double ik = sqrt(4*(z*p)/(1-a)); \
    ld = 2*(  \
            (1-5*z2+p2+q*q)*ellint_1<double>(ik) \
            + (1-a)*(z2+7*p2-4)*ellint_2<double>(ik) \
            - 3*q/a*ellint_3<double, double>(-ik, (a-b)/a) \
           ) / (9*M_PI*sqrt(1-a)); \
}

#define LAMBDA_3 { \
    double hip = 0.5/p; \
    ld = 1.0/3 + (16*p*(2*p2-1)*ellint_2<double>(hip) \
         - (1-4*p2)*(3-8*p2)/p*ellint_1<double>(hip) )/(9*M_PI); \
}

#define LAMBDA_4 { \
    double twop = 2*p; \
    ld = 1.0/3 + 2*( \
            4*(2*p2-1)*ellint_2<double>(twop) \
            + (1-4*p2)*ellint_1<double>(twop) \
         ) / (9*M_PI); \
}

#define LAMBDA_5 { \
    ld = 2*acos(1-2*p)/(3*M_PI) - 4*(3+2*p-8*p2)*sqrt(p*(1-p))/(9*M_PI); \
    if (z < p) ld -= 2.0/3; \
}

#define LAMBDA_6 { \
    ld = -2*(1-p2)*sqrt(1-p2)/3; \
}

#define ETA_2 { \
    eta = 0.5*p2*(p2+2*z2); \
}

#define ETA_1 { \
    ETA_2; \
    eta = (kap1+2*eta*kap0-0.25*(1+5*p2+z2)*sqrt((1-a)*(b-1)))/(2*M_PI); \
}

double min (double a, double b)
{
    if (a < b) return a;
    return b;
}

double max (double a, double b)
{
    if (a >= b) return a;
    return b;
}

double ldlc (double p, double z, double u1, double u2)
{
    double ld, le, eta, omega = 1.0 - u1/3.0 - u2/6.0,
           p2, z2, a, b, q, kap0, kap1, Kk, Ek, Pk;

    // Make sure that z is positive.
    z = fabs(z);

    if (p <= DBL_EPSILON || z >= 1+p) {
        ld = 0.0;
        eta = 0.0;
        le = 0.0;
    } else if (p >= 1.0 && z <= p-1.0) {
        ld = 0.0;
        eta = 0.5;
        le = 1.0;
    } else {
        // Pre-compute some common factors.
        a = (z-p)*(z-p);
        b = (z+p)*(z+p);
        p2 = p*p;
        z2 = z*z;
        q = p2 - z2;

        // Compute the non-limb-darkened transit.
        if (z >= fabs(1.0-p) && z <= 1.0+p) {
            kap1 = acos(min(0.5*(1-q)/z, 1));
            kap0 = acos(min(0.5*(p2+z2-1)/(p*z), 1));
            le = p2*kap0 + kap1 - 0.5*sqrt(max(4*z2-(1-q)*(1-q), 0.0));
            le /= M_PI;
        } else if (z <= 1-p) le = p2;

        // Compute the limb-darkened cases.
        if (0.5+fabs(p-0.5) < z && z < 1+p) {
#ifdef VERBOSE
            printf("Case 2\n");
#endif
            LAMBDA_1;
            ETA_1;
        } else if (p < 0.5 && p <= z && z <= 1-p) {
            if (z == 1-p) {
#ifdef VERBOSE
                printf("Case 4\n");
#endif
                LAMBDA_5;
            } else if (z == p) {
#ifdef VERBOSE
                printf("Case 5\n");
#endif
                LAMBDA_4;
            } else {
#ifdef VERBOSE
                printf("Case 3\n");
#endif
                LAMBDA_2;
            }
            ETA_2;
        } else if (p == 0.5 && z == 0.5) {
#ifdef VERBOSE
            printf("Case 6\n");
#endif
            ld = 1.0/3.0-4.0/(9*M_PI);
            eta = 3.0/32.0;
        } else if (p > 0.5 && fabs(1-p) <= z && z <= p) {
            if (z == p) {
#ifdef VERBOSE
                printf("Case 7\n");
#endif
                LAMBDA_3;
            } else if (z < p && z == 1-p) {
#ifdef VERBOSE
                printf("Case 12 - FIXME\n");
#endif
                // FIXME.
                ld = 0.0;
            } else {
#ifdef VERBOSE
                printf("Case 8\n");
#endif
                LAMBDA_1;
            }
            ETA_1;
        } else if (p < 1.0 && 0 <= z && z < 0.5-fabs(p-0.5)) {
            if (z < DBL_EPSILON) {
#ifdef VERBOSE
                printf("Case 10\n");
#endif
                LAMBDA_6;
            } else {
#ifdef VERBOSE
                printf("Case 9\n");
#endif
                LAMBDA_2;
            }
            ETA_2;
        } else {
            printf("FAIL\n");
            return 9999999999.0;
        }
    }

    if (z < p) ld += 2.0/3;

    return 1.0-((1-u1-2*u2)*le + (u1+2*u2)*ld + u2*eta)/omega;
}
