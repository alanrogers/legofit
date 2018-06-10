#ifndef ARR_HESSIAN_H
#define ARR_HESSIAN_H
#include <gsl/gsl_matrix.h>

typedef struct Hessian Hessian;

struct Hessian {
    int npar;
    char **parname;
    gsl_matrix *hessian;
};

Hessian hessian(const char *fname);

#endif
