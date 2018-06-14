#ifndef BRANCHTAB_INCLUDED
#  define BRANCHTAB_INCLUDED

#include "typedefs.h"
#include <stdio.h>

BranchTab    *BranchTab_new(void);
void          BranchTab_free(BranchTab * self);
double        BranchTab_get(BranchTab * self, tipId_t tipid);
int           BranchTab_hasSingletons(BranchTab * self);
void          BranchTab_add(BranchTab * self, tipId_t key, double value);
unsigned      BranchTab_size(BranchTab * self);
void          BranchTab_print(const BranchTab *self, FILE *fp);
void          BranchTab_plusEquals(BranchTab *lhs, BranchTab *rhs);
void          BranchTab_toArrays(BranchTab *self, unsigned n, tipId_t key[n],
                                 double value[n], double sqr[n]);
int           BranchTab_divideBy(BranchTab *self, double denom);
BranchTab    *BranchTab_parse(const char *fname, const LblNdx *lblndx);
BranchTab    *BranchTab_dup(const BranchTab *old);
int           BranchTab_equals(const BranchTab *lhs, const BranchTab *rhs);
double        BranchTab_sum(const BranchTab *self);
double        BranchTab_entropy(const BranchTab *self);
int           BranchTab_normalize(BranchTab *self);
double        BranchTab_chiSqCost(const BranchTab *obs, const BranchTab *expt,
                             double u, long nnuc, double n);
double        BranchTab_smplChiSqCost(const BranchTab *obs,
                                      const BranchTab *expt,
                                      double u, long nnuc, double n);
double        BranchTab_KLdiverg(const BranchTab *obs, const BranchTab *expt);
double        BranchTab_negLnL(const BranchTab *obs, const BranchTab *expt);
double        BranchTab_poissonCost(const BranchTab *obs, const BranchTab *expt,
                                    double u, long nnuc, double n);
#endif
