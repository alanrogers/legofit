#ifndef ARR_MCTREE_H
#  define ARR_MCTREE_H

#  include "typedefs.h"
#  include <stdio.h>
#  include <gsl/gsl_rng.h>

void       *MCTree_dup(const void *old);
int         MCTree_equals(const MCTree * lhs, const MCTree * rhs);
int         MCTree_feasible(const void * vself, int verbose);
void        MCTree_free(void *self);
LblNdx      MCTree_getLblNdx(void *self);
const char *MCTree_getNameFree(void * self, int i);
void        MCTree_getParams(void *self, int n, double x[n]);
void       *MCTree_new(const char *fname, Bounds bnd);
int         MCTree_nFree(const void *self);
void        MCTree_brlen(void * vself, BranchTab * branchtab, gsl_rng * rng,
                         unsigned long nreps, int doSing,
                         long unsigned *event_counter);
void        MCTree_printParStore(void * vself, FILE * fp);
void        MCTree_printParStoreFree(void *self, FILE *fp);
void        MCTree_randomize(void *self, gsl_rng *rng);
void        MCTree_sanityCheck(void * vself, const char *file, int line);
int         MCTree_setParams(void *self, int n, double x[n]);
void        MCTree_initStateVec(void *gpt, int ndx, int n, double x[n],
                                gsl_rng *rng);
const char *MCTree_getNameFree(void * self, int i);
void        MCTree_getParams(void * self, int n, double x[n]);
int         MCTree_nFree(const void * self);
int         MCTree_setParams(void * vself, int n, double x[n]);
void        MCTree_initStateVec(void *self, int ndx, int n, double x[n],
                                gsl_rng *rng);
unsigned    MCTree_nSamples(void *vself);
#endif
