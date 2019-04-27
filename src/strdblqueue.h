#ifndef STR_DBL_QUEUE
#define STR_DBL_QUEUE

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

typedef struct StrDblQueue StrDblQueue;
typedef struct StrDbl StrDbl;

struct StrDbl {
    char str[100];
    double val;
};

// A FIFO queue. New values are pushed onto the tail. Old ones are
// popped off of the head.
struct StrDblQueue {
    struct StrDblQueue *next;
    struct StrDbl strdbl;
};

StrDblQueue *StrDblQueue_push(StrDblQueue *prev, const char *str, double val);
StrDblQueue *StrDblQueue_pop(StrDblQueue *self, StrDbl *strdbl);
StrDblQueue *StrDblQueue_free(StrDblQueue *self);
int          StrDblQueue_length(StrDblQueue *self);
void         StrDblQueue_print(const StrDblQueue *self, FILE *fp);
int          StrDblQueue_compare(StrDblQueue *lhs, StrDblQueue *rhs);
// StrDbl      *StrDblQueue_get(StrDblQueue *self, StrDbl *strdbl, int index);
StrDblQueue *StrDblQueue_parseLegofit(const char *fname);
StrDblQueue *StrDblQueue_parseSitePat(const char *fname);
double StrDblQueue_msd(const StrDblQueue *a, const StrDblQueue *b);
void StrDblQueue_normalize(StrDblQueue *self);
void make_covar_matrix(int nfiles, int npar, double array[nfiles][npar],
                      gsl_matrix *cov);

#endif
