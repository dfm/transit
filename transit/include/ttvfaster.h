#ifndef _TRANSIT_TTVFASTER_H_
#define _TRANSIT_TTVFASTER_H_

//
// This is an adaptation of the TTVFaster C-implementation written by
// Katherine Deck available at https://github.com/ericagol/TTVFaster under the
// MIT license.
//

#include <cmath>
#include <vector>

#include "ceres/jet.h"

using std::abs;
using std::vector;
using ceres::Jet;

namespace transit {
namespace ttvfaster {

/*Many thanks for Jack Wisdom at MIT for sharing the code "laplace()" below with us. */
/* the longitude of transit = 0 */
/* Please cite Agol & Deck (2015) if you make use of this code in your research.*/
#define LAPLACE_EPS   1.0e-10
#define PI            M_PI
#define TWOPI         (2*M_PI)
#define MAX_N_PLANETS 9

// Note: - this defines the unit system: Msun, Day, AU
//       - this is different from the rest of the transit namespace.
#define G             2.95994511e-4;

// Gradient helpers
int ttv_to_int (double x) { return int(x); }
template <typename T, int N>
int ttv_to_int (Jet<T, N> x) { return int(x.a); }

double principal_value(double theta) {
    theta -= 2.0*PI*floor(theta/(2.0*PI));
    return(theta);
}
template <typename T, int N>
Jet<T, N> principal_value(Jet<T, N> theta) {
    double value = principal_value(theta.a);
    Jet<T, N> result(value, theta.v);
    return result;
}

// This is a wrapper around the results dataset.
template <typename T>
class TTVFasterResult {
public:
    TTVFasterResult (unsigned nplanets) : nplanets_(nplanets) {
        times_ = new T*[nplanets_];
        ntransits_ = new unsigned[nplanets_];
    };
    ~TTVFasterResult () {
        delete ntransits_;
        for (unsigned i = 0; i < nplanets_; ++i) delete times_[i];
        delete times_;
    };

    void set_ntransits (unsigned i, unsigned n) {
        times_[i] = new T[n];
        ntransits_[i] = n;
    };
    unsigned get_ntransits (unsigned i) const { return ntransits_[i]; };
    T* get_times (unsigned i) { return times_[i]; };

private:
    unsigned nplanets_, *ntransits_;
    T** times_;
};

template <typename T>
class TTVFaster {
public:

TTVFaster () {};
~TTVFaster () {};

// Constants for later use.
T AJ00, AJ01, AJ02, AJ10, AJ20, AJ11;

// Coefficient helper functions
T laplace(T s, int i, int j, T a) {
    T as, term, sum, factor1, factor2, factor3, factor4;
    int k,q, q0;

    as = a*a;

    if(i<0) i = -i;

    if(j<=i)     /* compute first term in sum */
    {
        factor4 = T(1.0);
        for(k=0; k<j; k++)
            factor4 = factor4 * double(i - k);
        sum = T(factor4);
        q0=0;
    }
    else
    {
        q0 = (j + 1 - i) / 2;
        sum = T(0.0);
        factor4 = T(1.0);
    }

    /* compute factors for terms in sum */

    factor1 = s;
    factor2 = s + double(i);
    factor3 = T(i + 1.0);
    for(q=1;q<q0;q++)   /* no contribution for q = 0 */
    {
        factor1 *= s + double(q);
        factor2 *= s + double(i + q);
        factor3 *= T(i + 1.0 + q);
    }

    term = as * factor1 * factor2 / (factor3 * double(q));

    /* sum series */

    while(term*factor4 > LAPLACE_EPS)
    {
        factor4 = T(1.0);
        for(k=0;k<j;k++)
            factor4 *= T(2.0*q + i - k);
        sum += term * factor4;
        factor1 += T(1.0);
        factor2 += T(1.0);
        factor3 += T(1.0);
        q++;
        term *= as * factor1 * factor2 / (factor3 * double(q));

    }

    /* fix coefficient */

    for(k=0;k<i;k++)
        sum *= (s + (double(k)))/((double(k))+1.0);

    if(q0 <= 0)
        sum *= 2.0 * pow(a, (double(i)));
    else
        sum *= 2.0 * pow(a, (double(2*q0 + i - 2)));

    return(sum);
}

T gam(int i, int k, int j, T kappa, T beta, T alpha) {
    T val = T(0.0);
    if(i==1){
        if(k ==0){
            val = beta;
        }
        if(k == 1){
            val = beta+1.0;
        }
        if(k == -1){
            val = beta-1.0;
        }
        if(k == 2){
            val = beta+pow(alpha,1.5);;
        }
        if(k == -2){
            val = beta-pow(alpha,1.5);
        }
    }else{
        if(k ==0){
            val = kappa;
        }
        if(k == 2){
            val = kappa+1.0;
        }
        if(k == -2){
            val = kappa-1.0;
        }
        if(k == 1){
            val = kappa+pow(alpha,-1.5);;
        }
        if(k == -1){
            val = kappa-pow(alpha,-1.5);
        }
    }
    return(val);
}

T xi(int i, int k, int j, T kappa, T beta, T alpha) {
    if(i==1) return beta;
    else return kappa;
}

T d1(int i, int k, int j, T kappa, T beta, T alpha) {
    T val = T(0.0);
    if(i==1){
        if(j==1){
            val = alpha*double(j)*(AJ00-alpha);
        }else{
            val = alpha*double(j)*AJ00;
        }
    }else{
        if(j==1){
            val = -double(j)*(AJ00-pow(alpha,-2.0));
        }else{
            val = -double(j)*AJ00;
        }
    }

    return(val);
}

T d2(int i, int k, int j, T kappa, T beta, T alpha) {
    T val = T(0.0);
    if(i==1){
        if(j==1){
            val = alpha*(AJ10-alpha);
        }else{
            val = alpha*AJ10;
        }
    }else{
        if(j==1){
            val = AJ01-pow(alpha,-2.0);
        }else{
            val = AJ01;
        }
    }

    return(val);
}

T C1(int i, int k, int j, T kappa, T beta, T alpha) {
    T val=T(0.0);
    if(i==1){
        if(k==0){
            val = d1(i,k,j,kappa,beta,alpha);
        }
        if(k ==1){
            if(j==1){
                val = alpha*double(j)*(double(j)*AJ00-0.5*AJ10+0.5*(1.0-2.0)*alpha);
            }else{
                val = alpha*double(j)*(double(j)*AJ00-0.5*AJ10);
            }
        }
        if(k ==-1){
            if(j==1){
                val = alpha*double(j)*(-double(j)*AJ00-0.5*AJ10+0.5*(1.0+2.0)*alpha);
            }else{
                val = alpha*double(j)*(-double(j)*AJ00-0.5*AJ10);
            }
        }


        if(k ==2){
            val = alpha*double(j)*(-double(j)*AJ00-0.5*AJ01);
        }

        if(k ==-2){
            if(j==1){
                val = alpha*double(j)*(double(j)*AJ00-0.5*AJ01-2.0*alpha);
            }else{
                val = alpha*double(j)*(double(j)*AJ00-0.5*AJ01);
            }
        }
    }else{
        /* Outer Planet */
        if(k==0){
            val = d1(i,k,j,kappa,beta,alpha);
        }
        if(k ==2){
            if(j==1){
                val = -double(j)*(-double(j)*AJ00-0.5*AJ01+0.5*3.0*pow(alpha,-2.0));
            }else{
                val = -double(j)*(-double(j)*AJ00-0.5*AJ01);
            }
        }
        if(k ==-2){
            if(j==1){
                val = -double(j)*(double(j)*AJ00-0.5*AJ01+0.5*(1.0-2.0)*pow(alpha,-2.0));
            }else{
                val = -double(j)*(double(j)*AJ00-0.5*AJ01);
            }
        }

        if(k ==-1){
            val = -double(j)*(-double(j)*AJ00-0.5*AJ10);
        }

        if(k ==1){
            if(j==1){
                val = -double(j)*(double(j)*AJ00-0.5*AJ10-2.0*pow(alpha,-2.0));

            }else{
                val = -double(j)*(double(j)*AJ00-0.5*AJ10);

            }
        }
    }

    return(val);
}

T C2(int i, int k, int j, T kappa, T beta, T alpha) {
    T val = T(0.0);
    if(i==1){
        if(k==0){
            val = d2(i,k,j,kappa,beta,alpha);
        }
        if(k ==1){
            if(j==1){
                val = alpha*(double(j)*AJ10-0.5*AJ20-alpha);
            }else{
                val = alpha*(double(j)*AJ10-0.5*AJ20);
            }
        }
        if(k ==-1){
            if(j==1){
                val = alpha*(-double(j)*AJ10-0.5*AJ20+alpha);
            }else{
                val = alpha*(-double(j)*AJ10-0.5*AJ20);
            }
        }

        if(k ==2){
            val = alpha*(-double(j)*AJ10-0.5*AJ11);
        }

        if(k ==-2){
            if(j==1){
                val = alpha*(double(j)*AJ10-0.5*AJ11-2.0*alpha);
            }else{
                val = alpha*(double(j)*AJ10-0.5*AJ11);
            }
        }
    }else{
        /* Outer Planet */
        if(k==0){
            val = d2(i,k,j,kappa,beta,alpha);
        }
        if(k ==2){
            if(j==1){
                val = (-double(j)*AJ01-0.5*AJ02+pow(alpha,-2.0));
            }else{
                val = (-double(j)*AJ01-0.5*AJ02);
            }
        }
        if(k ==-2){
            if(j==1){
                val = (double(j)*AJ01-0.5*AJ02-pow(alpha,-2.0));
            }else{
                val = (double(j)*AJ01-0.5*AJ02);
            }
        }

        if(k ==-1){
            val = (-double(j)*AJ01-0.5*AJ11);
        }

        if(k ==1){
            if(j==1){
                val = (double(j)*AJ01-0.5*AJ11-2.0*pow(alpha,-2.0));

            }else{
                val = (double(j)*AJ01-0.5*AJ11);

            }
        }
    }

    return(val);
}

T u(T gamma, T C1v, T C2v){
    T val=T(0.0);
    T gsq = gamma*gamma;
    val = ((3.0+gsq)*C1v+2.0*gamma*C2v)/gsq/(1.0-gsq);
    return(val);
}

T v_plus(T z, T d1v, T d2v){
    T val=T(0.0);
    T zsq = z*z;
    val = (((1.0-zsq)+6.0*z)*d1v+(2.0+zsq)*d2v)/(z*(1.0-zsq)*(z+1.0)*(z+2.0));
    return(val);
}

T v_minus(T z, T d1v, T d2v){
    T val=T(0.0);
    T zsq = z*z;
    val = ((-(1.0-zsq)+6.0*z)*d1v+(2.0+zsq)*d2v)/(z*(1.0-zsq)*(z-1.0)*(z-2.0));
    return(val);
}

T f_j_k_i(T alpha,int j,int k, int i, T kappa_o_m, T beta_o_m,vector<T>& a00,vector<T>& a01, vector<T>& a10, vector<T>& a20, vector<T>& a02, vector<T>& a11) {
    T val=T(0.0);
    T beta,kappa;

    AJ00= a00[j];
    AJ01= a01[j];
    AJ02= a02[j];
    AJ20= a20[j];
    AJ10= a10[j];
    AJ11= a11[j];
    beta = beta_o_m*double(j);
    kappa = kappa_o_m*double(j);
    val=u(gam(i,k,j,kappa,beta,alpha),C1(i,k,j,kappa,beta,alpha),C2(i,k,j,kappa,beta,alpha));
    if(i == abs(k)){
        if( k > 0){
            val+= v_plus(xi(i,k,j,kappa,beta,alpha),d1(i,k,j,kappa,beta,alpha),d2(i,k,j,kappa,beta,alpha));
        }
        else{
            val+= v_minus(xi(i,k,j,kappa,beta,alpha),d1(i,k,j,kappa,beta,alpha),d2(i,k,j,kappa,beta,alpha));
        }
    }
    return(val);
}

// This is the main entry point to the code.
TTVFasterResult<T>* compute_times (
    unsigned n_planets,
    const T* params,
    double t0,
    double tf,
    unsigned m_max
) {
    T mstar = exp(params[0]);

    // Pre-allocate some coefficient arrays.
    vector<T> b(m_max+2),
              db(m_max+2),
              d2b(m_max+2),
              AJ00_arr(m_max+2),
              AJ10_arr(m_max+2),
              AJ01_arr(m_max+2),
              AJ02_arr(m_max+2),
              AJ20_arr(m_max+2),
              AJ11_arr(m_max+2);

    // Compute the number of transits for each planet.
    vector<unsigned> n_transits(n_planets);
    TTVFasterResult<T>* times = new TTVFasterResult<T>(n_planets);
    unsigned j;
    for (j = 0; j < n_planets; j++) {
        T period = exp(params[j*7+2]);
        n_transits[j] = ttv_to_int((tf - t0) / period + 1.0);
        times->set_ntransits(j, n_transits[j]);
    }

    // "Keplerian" times
    unsigned count;
    for (j = 0; j < n_planets; j++) {
        T P1 = exp(params[j*7+2]),
          TT1 = params[j*7+7];
        for (count = 0; count < n_transits[j]; count++)
            times->get_times(j)[count] = t0 + P1*double(count) + TT1;
    }

    // Now add TTVs; adjacent pairs only
    for(j = 0; j < n_planets-1; j++){
        unsigned i = j+1;

        // First planet parameters.
        T P1 = exp(params[j*7+2]),
          TT1 = params[j*7+7],
          m1 = exp(params[j*7+1]),
          e1 = params[j*7+6]*params[7*j+6]+params[j*7+3]*params[7*j+3],
          ap1 = principal_value(atan2(params[j*7+6],params[j*7+3])),
          Omega1 = principal_value(params[j*7+5]),
          trueAnom1 = principal_value(0.0-ap1-Omega1), /* at transit*/
          EccAnom1 = 2.0*atan(sqrt(1.0-e1)/sqrt(1.0+e1)*tan(trueAnom1/2.0)),
          MT1= EccAnom1-e1*sin(EccAnom1),
          n1 = 2.0*PI/P1;

        T P2 = exp(params[i*7+2]),
          TT2 = params[i*7+7],
          m2 = exp(params[i*7+1]),
          e2 = params[i*7+6]*params[7*i+6]+params[i*7+3]*params[7*i+3],
          ap2 = principal_value(atan2(params[i*7+6],params[i*7+3])),
          Omega2 = principal_value(params[i*7+5]),
          trueAnom2 = principal_value(0.0-ap2-Omega2), /* at transit*/
          EccAnom2 =  2.0*atan(sqrt(1.0-e2)/sqrt(1.0+e2)*tan(trueAnom2/2.0)),
          MT2= EccAnom2-e2*sin(EccAnom2),
          n2 = 2.0*PI/P2;

        // ...
        T alpha = pow(P1/P2, 2.0/3.0);

        for (count = 0; count < m_max+2; count++){
            b[count] = laplace(T(0.5),count,0,alpha);
            db[count] = laplace(T(0.5),count,1,alpha);
            d2b[count] = laplace(T(0.5),count,2,alpha);
        }

        T kappa_over_m = (pow(alpha,-1.5)-1.0),
          beta_over_m = (1.0-pow(alpha,1.5));
        unsigned m;
        for (m = 0; m < m_max+2; m++) {
            AJ00_arr[m] = b[m];
            AJ10_arr[m] = db[m];
            AJ20_arr[m] = d2b[m];
            AJ01_arr[m] = -db[m]-b[m];
            AJ02_arr[m] = 4.0*db[m]+2.0*b[m]+d2b[m];
            AJ11_arr[m] = -2.0*db[m]-d2b[m];
        }

        /// Inner planet of pair
        for (count = 0; count < n_transits[j]; count++) {
            T base = P1*double(count) + TT1,
              lambda1 = n1*(base-TT1)+MT1+ap1+Omega1,
              lambda2 = n2*(base-TT2)+MT2+ap2+Omega2,
              psi = lambda1-lambda2,
              fac1 = lambda1-ap1-Omega1,
              fac2 = lambda1-ap2-Omega2,
              TTV = T(0.0);
            for (m = 1; m <= m_max; m++) {
                TTV += f_j_k_i(alpha,m,0,1,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(double(m)*psi);
                TTV += e1*f_j_k_i(alpha,m,1,1,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(double(m)*psi+fac1);
                TTV += e1*f_j_k_i(alpha,m,-1,1,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(double(m)*psi-fac1);
                TTV += e2*f_j_k_i(alpha,m+1,2,1,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(double(m)*psi+fac2);
                TTV += e2*f_j_k_i(alpha,m-1,-2,1,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(double(m)*psi-fac2);
            }

            TTV *= m2/mstar/n1;
            times->get_times(j)[count] += TTV;
        }

        // Outer planet of pair
        for (count = 0; count < n_transits[i]; count++) {
            T base = P2*double(count) + TT2,
              lambda1 = n1*(base-TT1)+MT1+ap1+Omega1,
              lambda2 = n2*(base-TT2)+MT2+ap2+Omega2,
              psi = lambda1-lambda2,
              fac1 = lambda2-ap1-Omega1,
              fac2 = lambda2-ap2-Omega2,
              TTV = T(0.0);

            for (m = 1; m <= m_max; m++) {
                TTV += f_j_k_i(alpha,m,0,2,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(double(m)*psi);
                TTV += e1*f_j_k_i(alpha,m-1,1,2,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(double(m)*psi+fac1);
                TTV += e1*f_j_k_i(alpha,m+1,-1,2,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(double(m)*psi-fac1);
                TTV += e2*f_j_k_i(alpha,m,2,2,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(double(m)*psi+fac2);
                TTV += e2*f_j_k_i(alpha,m,-2,2,kappa_over_m,beta_over_m,AJ00_arr,AJ01_arr,AJ10_arr,AJ20_arr,AJ02_arr,AJ11_arr)*sin(double(m)*psi-fac2);
            }

            TTV *= m1/mstar/n2;
            times->get_times(i)[count] += TTV;
        }
    }

    return times;
}

}; // TTVFaster

}; // namespace ttvfaster
}; // namespace transit

#endif // _TRANSIT_TTVFASTER_H_
