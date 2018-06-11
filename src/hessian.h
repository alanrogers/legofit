#ifndef ARR_HESSIAN_H
#define ARR_HESSIAN_H
#include <gsl/gsl_matrix.h>

typedef struct Hessian Hessian;

struct Hessian {
    int npar;            // number of parameters
    double lnL;          // log likelihood of estimate
    char **parname;      // array of pointers to parameter names
    gsl_matrix *hessian;
};

Hessian hessian(const char *fname);

#endif
