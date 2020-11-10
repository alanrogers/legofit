#ifndef ARR_GPTREE_H
#  define ARR_GPTREE_H

#  include "typedefs.h"
#  include <stdio.h>
#  include <gsl/gsl_rng.h>

void       *GPTree_dup(const void *old);
int         GPTree_feasible(const void * vself, int verbose);
void        GPTree_free(void * vself);
LblNdx      GPTree_getLblNdx(void *vself);
const char *GPTree_getNameFree(void * vself, int i);
void        GPTree_getParams(void *vself, int n, double x[n]);
void       *GPTree_new(const char *fname, Bounds bnd);
int         GPTree_nFree(const void *vself);
void        GPTree_brlen(void * vself, BranchTab * branchtab, gsl_rng * rng,
                         unsigned long nreps, int doSing, long unsigned *event_counter);
void        GPTree_patprob(void *vself, BranchTab *branchtab,
                            gsl_rng *rng, unsigned long nreps,
                            int doSing);
void        GPTree_printParStore(void *vself, FILE *fp);
void        GPTree_printParStoreFree(void *vself, FILE *fp);
void        GPTree_randomize(void *vself, gsl_rng *rng);
void        GPTree_sanityCheck(void *vself, const char *file, int line);
int         GPTree_setParams(void *vself, int n, double x[n]);
void        GPTree_initStateVec(void *gpt, int ndx, int n, double x[n],
                                gsl_rng *rng);
unsigned    GPTree_nSamples(void *vself);
#endif


