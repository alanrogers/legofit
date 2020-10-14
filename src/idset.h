#ifndef ARR_IDSET_H
#  define ARR_IDSET_H

#  include "typedefs.h"
#  include "eventlst.h"
#  include <stdio.h>
#  include <stdlib.h>

IdSet *IdSet_addSamples(IdSet *old, int nsamples, tipId_t *samples);
void   IdSet_addEvent(IdSet *self, unsigned event, unsigned outcome,
                      long double pr);
IdSet *IdSet_join(IdSet *left, IdSet *right, int nsamples,
                  tipId_t *samples);
void   IdSet_copyEventLst(IdSet *self, const IdSet *old);
IdSet *IdSet_dup(const IdSet *old);
EventLst *IdSet_dupEventLst(const IdSet *self);
IdSet *IdSet_join(IdSet *left, IdSet *right, int nsamples,
                  tipId_t samples[nsamples]);
IdSet *IdSet_new(int nIds, const tipId_t *tid, EventLst *evlst);
IdSet *IdSet_newTip(tipId_t tid);
void   IdSet_print(IdSet *self, FILE *fp);
void   IdSet_sanityCheck(IdSet *self, const char *file, int lineno);
void   IdSet_free(IdSet *self);
int    IdSet_nIds(IdSet *self);
long double IdSet_prob(IdSet *self);
int    IdSet_cmp(const IdSet *left, const IdSet *right);
uint32_t IdSet_hash(const IdSet *self) __attribute__((no_sanitize("integer")));
tipId_t IdSet_union(IdSet *self, int n, int *ndx);
void   IdSet_partition(IdSet *self, unsigned k, tipId_t part[k],
                       unsigned n, unsigned a[n]);
IdSet *IdSet_subset(const IdSet *self, int n, int *ndx);
int    IdSet_get(const IdSet *self, int i);

#endif
