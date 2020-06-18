/**
 * @file idset.c
 * @author Alan R. Rogers
 * @brief A set of tipID_t values.
 *
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "idset.h"
#include "misc.h"
#include "binary.h"
#include <stdio.h>
#include <stdlib.h>

void IdSet_sanityCheck(IdSet *self, const char *file, int lineno) {
#ifndef NDEBUG
    if(self == NULL)
       return;
    REQUIRE(self->p >= 0.0, file, lineno);
    REQUIRE(self->p <= 1.0, file, lineno);
#endif    
}

void IdSet_print(IdSet *self, FILE *fp) {
    if(self==NULL) {
        fputs("-----------------------\n", fp);
        return;
    }
    fprintf(fp, "probability = %lf nIds=%d\n",
            self->p, self->nIds);
    for(int i=0; i < self->nIds; ++i)
        fprintf(fp, " 0%o", self->tid[i]);
    putc('\n', fp);
    fprintf(fp,"mig history:");
    MigOutcome_print(self->mig, fp);
}

/**
 * Compare IdSet objects. Return -1 if x<y, 1 if x>y, and 0 otherwise.
 */
int IdSet_cmp(const IdSet *x, const IdSet *y) {
    for(int i=0; i < x->nIds && i < y->nIds; ++i) {
        if(x->tid[i] < y->tid[i])
            return -1;
        if(x->tid[i] > y->tid[i])
            return 1;
    }
    if(x->nIds < y->nIds)
        return -1;
    if(x->nIds > y->nIds)
        return 1;
    return 0;
}

/// Allocate a new IdSet with a single tipId_t value and probability 1.
IdSet *IdSet_newTip(tipId_t tid) {
    IdSet *self = malloc(sizeof(IdSet));
    CHECKMEM(self);

    self->nIds = 1;
    self->p = 1.0;
    self->mig = NULL;
    self->tid[0] = tid;
    return self;
}

/**
 * Allocate a new IdSet object with given values of nIds, tid, and
 * prob. Uses the struct hack.
 */
IdSet *IdSet_new(int nIds, const tipId_t tid[nIds], double prob) {
    assert(nIds > 0);
    size_t size = sizeof(IdSet) + (nIds-1) * sizeof(tipId_t);
    
    IdSet *self = malloc(size);
    CHECKMEM(self);

    self->nIds = nIds;
    self->p = prob;
    self->mig = NULL;

    for(int i=0; i < nIds; ++i)
        self->tid[i] = tid[i];

    return self;
}

// Duplicate
IdSet *IdSet_dup(const IdSet *old) {
    if(old == NULL)
        return NULL;
    IdSet *new = IdSet_new(old->nIds, old->tid, old->p);
    new->mig = MigOutcome_dup(old->mig);
    return new;
}

// Add MigOutcome list from old into self.
void IdSet_copyMigOutcome(IdSet *self, const IdSet *old) {
    assert(self->mig == NULL);
    self->mig = MigOutcome_dup(old->mig);
}
