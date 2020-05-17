#ifndef ARR_MCTREE_H
#  define ARR_MCTREE_H

#  include "typedefs.h"
#  include <stdio.h>
#  include <gsl/gsl_rng.h>

MCTree     *MCTree_dup(const MCTree *old);
int         MCTree_feasible(const MCTree *self, int verbose);
void        MCTree_free(MCTree *self);
LblNdx      MCTree_getLblNdx(MCTree *self);
const char *MCTree_getNameFree(MCTree * self, int i);
void        MCTree_getParams(MCTree *self, int n, double x[n]);
MCTree     *MCTree_new(const char *fname, Bounds bnd);
int         MCTree_nFree(const MCTree *self);
void        MCTree_patprob(MCTree *self, BranchTab *branchtab,
                            gsl_rng *rng, unsigned long nreps,
                            int doSing);
void        MCTree_printParStore(MCTree *self, FILE *fp);
void        MCTree_printParStoreFree(MCTree *self, FILE *fp);
void        MCTree_randomize(MCTree *self, gsl_rng *rng);
void        MCTree_sanityCheck(MCTree *self, const char *file, int line);
int         MCTree_setParams(MCTree *self, int n, double x[n]);
void        MCTree_initStateVec(MCTree *gpt, int ndx, int n, double x[n],
                                gsl_rng *rng);
#endif


