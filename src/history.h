#ifndef ARR_HISTORY_H
#  define ARR_HISTORY_H

#  include "typedefs.h"
#  include "parstore.h"
#  include <stdio.h>
#  include <gsl/gsl_rng.h>

History    *History_dup(const History *old);
int         History_feasible(const History *self, int verbose);
void        History_free(History *self);
LblNdx      History_getLblNdx(History *self);
const char *History_getNameFree(History * self, int i);
void        History_getParams(History *self, int n, double x[n]);
History    *History_new(const char *fname, Bounds bnd, int method);
int         History_nFree(const History *self);
void        History_patprob(History *self, BranchTab *branchtab,
                            gsl_rng *rng, unsigned long nreps,
                            int doSing);
void        History_printParStore(History *self, FILE *fp);
void        History_printParStoreFree(History *self, FILE *fp);
void        History_sanityCheck(History *self, const char *file, int line);
int         History_setParams(History *self, int n, double x[n]);
void        History_initStateVec(History *self, int ndx, int n, double x[n],
                                 gsl_rng *rng);
void       *History_newNode(History *self, double *twoN, double *start,
                            NodeStore *ns);
int         History_addChild(History *self, void * parent, void * child);

#endif


