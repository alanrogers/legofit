#ifndef ARR_IDSET_H
#  define ARR_IDSET_H

#  include "typedefs.h"
#  include <stdio.h>

IdSet *IdSet_newTip(tipId_t tid);
IdSet *IdSet_new(IdSet *next, int nIds, tipId_t tid[nIds], double prob);
int    IdSet_cmp(const IdSet *x, const IdSet *y);
IdSet *IdSet_dup(const IdSet *old);
IdSet *IdSet_add(IdSet *head, const IdSet *to_add, double prob);
void   IdSet_free(IdSet *self);
IdSet *IdSet_join(IdSet *a, IdSet *b);
int    IdSet_nSets(IdSet *self);
int    IdSet_nDescendants(IdSet *self);
void   IdSet_mulBy(IdSet *self, double factor);
void   IdSet_divBy(IdSet *self, double divisor);
double IdSet_sumProb(IdSet *self);
void   IdSet_print(IdSet *self, FILE *fp);

#endif
