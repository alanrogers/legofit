#ifndef ARR_GPTREE_H
#  define ARR_GPTREE_H

#  include "typedefs.h"
#  include <stdio.h>
#  include <gsl/gsl_rng.h>

GPTree     *GPTree_new(const char *fname, Bounds bnd);
void        GPTree_free(GPTree *self);
void        GPTree_sanityCheck(GPTree *self, const char *file, int line);
LblNdx      GPTree_getLblNdx(GPTree *self);
void        GPTree_patprob(GPTree *self, BranchTab *branchtab,
                            gsl_rng *rng, unsigned long nreps,
                            int doSing);
int         GPTree_nFree(const GPTree *self);
int         GPTree_setParams(GPTree *self, int n, double x[n]);
void        GPTree_getParams(GPTree *self, int n, double x[n]);
void        GPTree_printParStore(GPTree *self, FILE *fp);
void        GPTree_printParStoreFree(GPTree *self, FILE *fp);
const char *GPTree_getNameFree(GPTree * self, int i);
int         GPTree_feasible(const GPTree *self, int verbose);
void        initStateVec(int ndx, GPTree *gpt, int n, double x[n],
                         gsl_rng *rng);

// Not used outside of gptree.c and popnode.c
GPTree     *GPTree_dup(const GPTree *old);
int         GPTree_equals(const GPTree *lhs, const GPTree *rhs);
unsigned    GPTree_nsamples(GPTree *self);
void        GPTree_randomize(GPTree *self, gsl_rng *rng);
#endif


