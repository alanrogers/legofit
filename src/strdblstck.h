#ifndef STR_DBL_STACK
#define STR_DBL_STACK

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

typedef struct StrDblStack StrDblStack;
typedef struct StrDbl StrDbl;

struct StrDbl {
    char str[100];
    double val;
};

// A FIFO stack. New values are pushed onto the tail. Old ones are
// popped off of the head.
struct StrDblStack {
    struct StrDblStack *next;
    struct StrDbl strdbl;
};

StrDblStack *StrDblStack_push(StrDblStack *prev, const char *str, double val);
StrDblStack *StrDblStack_pop(StrDblStack *self, StrDbl *strdbl);
StrDblStack *StrDblStack_free(StrDblStack *self);
int          StrDblStack_length(StrDblStack *self);
void         StrDblStack_print(StrDblStack *self, FILE *fp);
int          StrDblStack_compare(StrDblStack *lhs, StrDblStack *rhs);
StrDbl      *StrDblStack_get(StrDblStack *self, StrDbl *strdbl, int index);
StrDblStack *parseLegofit_CLIC(const char *fname);
StrDblStack *parseLegofit_BEPE(const char *fname);
StrDblStack *parseData_BEPE(const char *fname);
StrDblStack *normalize(StrDblStack *self);
void make_covar_matrix(int nfiles, int npar, double array[nfiles][npar],
                      gsl_matrix *cov);

#endif
