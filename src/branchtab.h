#ifndef BRANCHTAB_INCLUDED
#  define BRANCHTAB_INCLUDED

#include "typedefs.h"
#include <stdio.h>

BranchTab    *BranchTab_new(unsigned nsamples);
void          BranchTab_free(BranchTab * self);
long double   BranchTab_get(BranchTab * self, tipId_t tipid);
int           BranchTab_hasSingletons(BranchTab * self);
void          BranchTab_add(BranchTab * self, tipId_t key, long double value);
unsigned      BranchTab_size(BranchTab * self);
void          BranchTab_print(const BranchTab *self, FILE *fp);

void          BranchTab_plusEquals(BranchTab *lhs, BranchTab *rhs);
void          BranchTab_minusEquals(BranchTab *lhs, BranchTab *rhs);
void          BranchTab_toArrays(BranchTab *self, unsigned n, tipId_t key[n],
                                 long double value[n]);
int           BranchTab_divideBy(BranchTab *self, long double denom);
BranchTab    *BranchTab_dup(const BranchTab *old);
int           BranchTab_equals(const BranchTab *lhs, const BranchTab *rhs);
long double   BranchTab_sum(const BranchTab *self);
long double   BranchTab_entropy(const BranchTab *self);
int           BranchTab_normalize(BranchTab *self);
long double   BranchTab_chiSqCost(const BranchTab *obs, const BranchTab *expt,
                             long double u, long nnuc, long double n);
long double   BranchTab_smplChiSqCost(const BranchTab *obs,
                                      const BranchTab *expt,
                                      long double u, long nnuc, long double n);
long double   BranchTab_KLdiverg(const BranchTab *obs, const BranchTab *expt);
long double   BranchTab_negLnL(const BranchTab *obs, const BranchTab *expt);
long double   BranchTab_poissonCost(const BranchTab *obs, const BranchTab *expt,
                                    long double u, long nnuc, long double n);
BranchTab    *BranchTab_collapse(BranchTab *old, tipId_t collapse);
BranchTab    *BranchTab_rmPops(BranchTab *old, tipId_t remove);
void          BranchTab_sanityCheck(BranchTab *self, const char *file,
                                    int line);
void          make_collapse_map(size_t n, tipId_t map[n], tipId_t collapse);
void          make_rm_map(size_t n, tipId_t map[n], tipId_t remove);
tipId_t       remap_bits(size_t n, tipId_t map[n], tipId_t old);
#endif
