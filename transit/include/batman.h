#ifndef _TRANSIT_BATMAN_H_
#define _TRANSIT_BATMAN_H_

#include <cmath>
#include <cfloat>
#include "ellint.h"

using std::abs;

namespace transit {

#define MIN(x, y) (((x) < (y)) ? T(x) : T(y))
#define MAX(x, y) (((x) > (y)) ? T(x) : T(y))

class BatmanLimbDarkening {
public:
    BatmanLimbDarkening () {};

    template <typename T>
    T operator () (const T* const params, const T& p, const T& z) const
    {
        double tol = DBL_EPSILON;
        T q1 = sqrt(params[0]), q2 = 2.0 * params[1],
          c1 = q1 * q2, c2 = q1 * (1.0 - q2),
          omega = 1.0 - c1/3.0 - c2/6.0,
          lambdad = T(0.0), etad = T(0.0), lambdae = T(0.0),
          kap1 = T(0.0), kap0 = T(0.0),
          Kk, Ek, Pk, q, n, x1, x2, x3, d;

        // allow for negative impact parameters
        d = abs(z);

        // check the corner cases
        if(abs(p - d) < tol)
        {
            d = p;
        }
        if(abs(p - 1.0 - d) < tol)
        {
            d = p - 1.0;
        }
        if(abs(1.0 - p - d) < tol)
        {
            d = 1.0 - p;
        }
        if(d < tol)
        {
            d = T(0.0);
        }

        x1 = (p - d) * (p - d);
        x2 = (p + d) * (p + d);
        x3 = p*p - d*d;

        //source is unocculted:
        if(d >= 1.0 + p)
        {
            //printf("zone 1\n");
            return T(1.0);
        }
        //source is completely occulted:
        if(p >= 1.0 && d <= p - 1.0)
        {
            //printf("zone 2\n");
            lambdad = T(0.0);
            etad = T(0.5);        //error in Fortran code corrected here, following Jason Eastman's python code
            lambdae = T(1.0);
            return 1.0 - ((1.0 - c1 - 2.0*c2)*lambdae + (c1 + 2.0*c2)*(lambdad + 2.0/3.0) + c2*etad)/omega;
        }
        //source is partly occulted and occulting object crosses the limb:
        if(d >= abs(1.0 - p) && d <= 1.0 + p)
        {
            //printf("zone 3\n");
            kap1 = acos(MIN((1.0 - p*p + d*d)/2.0/d, 1.0));
            kap0 = acos(MIN((p*p + d*d - 1.0)/2.0/p/d, 1.0));
            lambdae = p*p*kap0 + kap1;
            lambdae = (lambdae - 0.50*sqrt(MAX(4.0*d*d - pow((1.0 + d*d - p*p), 2.0), 0.0)))/M_PI;
        }

        //edge of the occulting star lies at the origin
        if(d == p)
        {
            //printf("zone 5\n");
            if(d < 0.5)
            {
                //printf("zone 5.2\n");
                q = 2.0*p;
                Kk = ellint_2(q);
                Ek = ellint_1(q);
                lambdad = 1.0/3.0 + 2.0/9.0/M_PI*(4.0*(2.0*p*p - 1.0)*Ek + (1.0 - 4.0*p*p)*Kk);
                etad = p*p/2.0*(p*p + 2.0*d*d);
                return 1.0 - ((1.0 - c1 - 2.0*c2)*lambdae + (c1 + 2.0*c2)*lambdad + c2*etad)/omega;
            }
            else if(d > 0.5)
            {
                //printf("zone 5.1\n");
                q = 0.5/p;
                Kk = ellint_2(q);
                Ek = ellint_1(q);
                lambdad = 1.0/3.0 + 16.0*p/9.0/M_PI*(2.0*p*p - 1.0)*Ek -  \
                          (32.0*pow(p, 4.0) - 20.0*p*p + 3.0)/9.0/M_PI/p*Kk;
                etad = 1.0/2.0/M_PI*(kap1 + p*p*(p*p + 2.0*d*d)*kap0 -  \
                                  (1.0 + 5.0*p*p + d*d)/4.0*sqrt((1.0 - x1)*(x2 - 1.0)));
            //    continue;
            }
            else
            {
                //printf("zone 6\n");
                lambdad = T(1.0/3.0 - 4.0/M_PI/9.0);
                etad = T(3.0/32.0);
                return 1.0 - ((1.0 - c1 - 2.0*c2)*lambdae + (c1 + 2.0*c2)*lambdad + c2*etad)/omega;
            }

            return 1.0 - ((1.0 - c1 - 2.0*c2)*lambdae + (c1 + 2.0*c2)*lambdad + c2*etad)/omega;
        }
         //occulting star partly occults the source and crosses the limb:
        //if((d > 0.5 + abs(p  - 0.5) && d < 1.0 + p) || (p > 0.5 && d > abs(1.0 - p)*1.0001 \
        //&& d < p))  //the factor of 1.0001 is from the Mandel/Agol Fortran routine, but gave bad output for d near abs(1-p)
        if((d > 0.5 + abs(p  - 0.5) && d < 1.0 + p) || (p > 0.5 && d > abs(1.0 - p) \
            && d < p))
        {
            //printf("zone 3.1\n");
            q = sqrt((1.0 - x1)/4.0/d/p);
            Kk = ellint_2(q);
            Ek = ellint_1(q);
            n = 1.0/x1 - 1.0;
            Pk = ellint_3(-q, n);
            lambdad = 1.0/9.0/M_PI/sqrt(p*d)*(((1.0 - x2)*(2.0*x2 +  \
                    x1 - 3.0) - 3.0*x3*(x2 - 2.0))*Kk + 4.0*p*d*(d*d +  \
                    7.0*p*p - 4.0)*Ek - 3.0*x3/x1*Pk);
            if(d < p) lambdad += T(2.0/3.0);
            etad = 1.0/2.0/M_PI*(kap1 + p*p*(p*p + 2.0*d*d)*kap0 -  \
                (1.0 + 5.0*p*p + d*d)/4.0*sqrt((1.0 - x1)*(x2 - 1.0)));
            return 1.0 - ((1.0 - c1 - 2.0*c2)*lambdae + (c1 + 2.0*c2)*lambdad + c2*etad)/omega;
        }
        //occulting star transits the source:
        if(p <= 1.0  && d <= (1.0 - p))
        {
            etad = p*p/2.0*(p*p + 2.0*d*d);
            lambdae = p*p;

            //printf("zone 4.1\n");
            q = sqrt((x2 - x1)/(1.0 - x1));
            Kk = ellint_2(q);
            Ek = ellint_1(q);
            n = x2/x1 - 1.0;
            Pk = ellint_3(-q, n);

            lambdad = 2.0/9.0/M_PI/sqrt(1.0 - x1)*((1.0 - 5.0*d*d + p*p +  \
                     x3*x3)*Kk + (1.0 - x1)*(d*d + 7.0*p*p - 4.0)*Ek - 3.0*x3/x1*Pk);

            // edge of planet hits edge of star
            if(abs(p + d - 1.0) <= tol)
            {
                lambdad = 2.0/3.0/M_PI*acos(1.0 - 2.0*p) - 4.0/9.0/M_PI* \
                            sqrt(p*(1.0 - p))*(3.0 + 2.0*p - 8.0*p*p);
            }
            if(d < p) lambdad += T(2.0/3.0);
        }
        return 1.0 - ((1.0 - c1 - 2.0*c2)*lambdae + (c1 + 2.0*c2)*lambdad + c2*etad)/omega;
    };
};

};

#endif
