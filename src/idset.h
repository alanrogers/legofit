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

void   IdSet_addMigEvent(IdSet *self, unsigned event, unsigned outcome,
                         double pr);
IdSet *IdSet_join(IdSet *left, IdSet *right, int nsamples,
                  tipId_t *samples);
void   IdSet_copyMigOutcome(IdSet *self, const IdSet *old);
IdSet *IdSet_dup(const IdSet *old);
IdSet *IdSet_join(IdSet *left, IdSet *right, int nsamples,
                  tipId_t samples[nsamples]);
IdSet *IdSet_new(int nIds, const tipId_t tid[nIds], double prob);
IdSet *IdSet_newTip(tipId_t tid);
void   IdSet_print(IdSet *self, FILE *fp);
void   IdSet_sanityCheck(IdSet *self, const char *file, int lineno);
static inline void IdSet_free(IdSet *self);
static inline int IdSet_nIds(IdSet *self);

static inline void IdSet_free(IdSet *self) {
    MigOutcome_free(self->mig);
    free(self);
}
static inline int IdSet_nIds(IdSet *self) { return self->nIds; }

#endif
