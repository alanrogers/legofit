#ifndef ARR_IDSET_H
#  define ARR_IDSET_H

#  include "typedefs.h"
#  include <stdio.h>
#  include <stdlib.h>

// A set of tipId_t values.
struct IdSet {
    double p; // probability of this set
    int nIds; // number of Ids in this set

    // Array of length nIds. tid[i] is the id of the i'th haploid
    // individual in this set. Allocated using the struct hack.
    tipId_t tid[1];
};

IdSet *IdSet_newTip(tipId_t tid);
IdSet *IdSet_new(int nIds, const tipId_t tid[nIds], double prob);
int    IdSet_cmp(const IdSet *x, const IdSet *y);
void   IdSet_print(IdSet *self, FILE *fp);
void   IdSet_sanityCheck(IdSet *self, const char *file, int lineno);
static inline IdSet *IdSet_dup(const IdSet *old);
static inline void IdSet_free(IdSet *self);
static inline int IdSet_nIds(IdSet *self);

static inline void IdSet_free(IdSet *self) {
    free(self);
}

// Duplicate
static inline IdSet *IdSet_dup(const IdSet *old) {
    if(old == NULL)
        return NULL;
    return IdSet_new(old->nIds, old->tid, old->p);
}

static inline int IdSet_nIds(IdSet *self) {
    return self->nIds;
}

#endif
