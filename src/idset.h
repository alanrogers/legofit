#ifndef ARR_IDSET_H
#  define ARR_IDSET_H

#  include "typedefs.h"
#  include "event.h"
#  include <stdio.h>
#  include <stdlib.h>

// A set of tipId_t values.
struct IdSet {
    long double p; // probability of this set
    int nIds; // number of Ids in this set

    Event *evlst; // linked list of events, outcomes and probs.

    // Array of length nIds. tid[i] is the id of the i'th haploid
    // individual in this set. Allocated using the struct hack.
    tipId_t tid[1];
};

IdSet *IdSet_addSamples(IdSet *old, int nsamples, tipId_t *samples);
void   IdSet_addEvent(IdSet *self, unsigned event, unsigned outcome,
                         long double pr);
IdSet *IdSet_join(IdSet *left, IdSet *right, int nsamples,
                  tipId_t *samples);
void   IdSet_copyEvent(IdSet *self, const IdSet *old);
IdSet *IdSet_dup(const IdSet *old);
IdSet *IdSet_join(IdSet *left, IdSet *right, int nsamples,
                  tipId_t samples[nsamples]);
IdSet *IdSet_new(int nIds, const tipId_t *tid, long double prob);
IdSet *IdSet_newTip(tipId_t tid);
void   IdSet_print(IdSet *self, FILE *fp);
void   IdSet_sanityCheck(IdSet *self, const char *file, int lineno);
static inline void IdSet_free(IdSet *self);
static inline int IdSet_nIds(IdSet *self);
static inline long double IdSet_prob(IdSet *self);
int    IdSet_cmp(const IdSet *left, const IdSet *right);
uint32_t IdSet_hash(const IdSet *self) __attribute__((no_sanitize("integer")));

static inline void IdSet_free(IdSet *self) {
    Event_free(self->evlst);
    free(self);
}

static inline int IdSet_nIds(IdSet *self) { return self->nIds; }

static inline long double IdSet_prob(IdSet *self) {
    return self->p * Event_prob(self->evlst);
}


#endif
