#ifndef ARR_GPTREE_H
#  define ARR_GPTREE_H

#  include "typedefs.h"

GPTree     *GPTree_new(const char *fname, Bounds bnd);
void        GPTree_free(GPTree *self);
GPTree     *GPTree_dup(const GPTree *old);
void        GPTree_sanityCheck(GPTree *self, const char *file, int line);
int         GPTree_equals(GPTree *lhs, GPTree *rhs);
LblNdx     *GPTree_getLblNdxPtr(GPTree *self);
void        GPTree_simulate(GPTree *self, BranchTab *branchtab,
                            gsl_rng *rng, unsigned long nreps);
#endif
