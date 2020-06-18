#ifndef ARR_IDSET_H
#  define ARR_IDSET_H

#  include "typedefs.h"
#  include "migoutcome.h"
#  include <stdio.h>
#  include <stdlib.h>

// A set of tipId_t values.
struct IdSet {
    double p; // probability of this set
    int nIds; // number of Ids in this set

    MigOutcome *mig;

    // Array of length nIds. tid[i] is the id of the i'th haploid
    // individual in this set. Allocated using the struct hack.
    tipId_t tid[1];
};

IdSet *IdSet_newTip(tipId_t tid);
IdSet *IdSet_new(int nIds, const tipId_t tid[nIds], double prob);
int    IdSet_cmp(const IdSet *x, const IdSet *y);
void   IdSet_print(IdSet *self, FILE *fp);
void   IdSet_sanityCheck(IdSet *self, const char *file, int lineno);
IdSet *IdSet_dup(const IdSet *old);
void   IdSet_copyMigOutcome(IdSet *self, const IdSet *old);
static inline void IdSet_free(IdSet *self);
static inline int IdSet_nIds(IdSet *self);

static inline void IdSet_free(IdSet *self) {
    MigOutcome_free(self->mig);
    free(self);
}
static inline int IdSet_nIds(IdSet *self) { return self->nIds; }

#endif
