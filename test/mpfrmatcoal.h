#ifndef MPFRMATCOAL_INCLUDED
#  define MPFRMATCOAL_INCLUDED

#  include "typedefs.h"
#  include <stdio.h>
#  include <mpfr.h>

typedef struct MpfrVec MpfrVec;
typedef struct MpfrMatCoal MpfrMatCoal;

struct MpfrVec {
    unsigned    dim;
    mpfr_t     *x;
};

MpfrVec    *MpfrVec_new(unsigned dim, long double x[dim]);
void        MpfrVec_free(MpfrVec * self);
void        MpfrVec_get(MpfrVec *self, unsigned dim, long double x[dim]);
void        MpfrVec_set(MpfrVec *self, unsigned dim,
                        const long double x[dim]);
void        MpfrVec_print(MpfrVec *self, FILE *ofp);

MpfrMatCoal *MpfrMatCoal_new(unsigned nSamples);
void        MpfrMatCoal_free(MpfrMatCoal *self);
void        MpfrMatCoal_print(MpfrMatCoal * self);
void        MpfrMatCoal_project(MpfrMatCoal *self, int dim,
                                long double ans[dim], long double v);
void        MpfrMatCoal_ciLen(MpfrMatCoal *self, unsigned dim,
                              long double m[dim], long double dv);

void UTmatXvec(unsigned dim, mpfr_t *y, mpfr_t *A, unsigned *offset, mpfr_t *x);

#endif
