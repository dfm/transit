#ifndef _QUAD_H_
#define _QUAD_H_

#ifdef __cplusplus
extern "C" {
#endif

double ldlc (double p, double z, double u1, double u2);

#ifdef __cplusplus
}
#endif

namespace transit {
class QuadraticLimbDarkening {

public:

    QuadraticLimbDarkening (double u1, double u2)
        : u1_(u1), u2_(u2) {};

    double operator () (double p, double z) const;

private:

    double u1_, u2_;

};
}

#endif
