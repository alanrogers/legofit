#ifndef ARR_GPTREE_H
#  define ARR_GPTREE_H

#  include "typedefs.h"
#  include <stdio.h>
#  include <gsl/gsl_rng.h>

GPTree     *GPTree_dup(const GPTree *old);
int         GPTree_feasible(const GPTree *self, int verbose);
void        GPTree_free(GPTree *self);
LblNdx      GPTree_getLblNdx(GPTree *self);
const char *GPTree_getNameFree(GPTree * self, int i);
void        GPTree_getParams(GPTree *self, int n, double x[n]);
GPTree     *GPTree_new(const char *fname, Bounds bnd);
int         GPTree_nFree(const GPTree *self);
void        GPTree_patprob(GPTree *self, BranchTab *branchtab,
                            gsl_rng *rng, unsigned long nreps,
                            int doSing);
void        GPTree_printParStore(GPTree *self, FILE *fp);
void        GPTree_printParStoreFree(GPTree *self, FILE *fp);
void        GPTree_sanityCheck(GPTree *self, const char *file, int line);
int         GPTree_setParams(GPTree *self, int n, double x[n]);
void        initStateVec(int ndx, GPTree *gpt, int n, double x[n],
                         gsl_rng *rng);
#endif


