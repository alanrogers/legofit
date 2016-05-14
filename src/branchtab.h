#ifndef BRANCHTAB_INCLUDED
#  define BRANCHTAB_INCLUDED

#include "typedefs.h"

BranchTab    *BranchTab_new(void);
void          BranchTab_free(BranchTab * self);
double        BranchTab_get(BranchTab * self, tipId_t tipid);
void          BranchTab_add(BranchTab * self, tipId_t key, double value);
unsigned      BranchTab_size(BranchTab * self);
void          BranchTab_print(BranchTab *self);
void          BranchTab_plusEquals(BranchTab *lhs, BranchTab *rhs);
void          BranchTab_toArrays(BranchTab *self, unsigned n, tipId_t key[n],
                                 double value[n]);
double        BranchTab_sum(const BranchTab *self);
int           BranchTab_normalize(BranchTab *self);
BranchTab    *BranchTab_parse(const char *fname, const LblNdx *lblndx);
BranchTab    *BranchTab_dup(const BranchTab *old);
int           BranchTab_equals(const BranchTab *lhs, const BranchTab *rhs);
double        BranchTab_KLdiverg(const BranchTab *obs, const BranchTab *expt);
#endif
