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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void IdSet_sanityCheck(IdSet *self, const char *file, int lineno) {
#ifndef NDEBUG
    if(self == NULL)
       return;
    REQUIRE(self->nIds >- 0, file, lineno);
    REQUIRE(self->p >= 0.0, file, lineno);
    REQUIRE(self->p <= 1.0, file, lineno);
    for(int i=0; i < self->nIds; ++i)
        REQUIRE(self->tid[i] > 0, file, lineno);

    // The bits set in each tipId_t value should be mutually
    // exclusive.
    tipId_t u = 0; // union of prior tipId_t values
    for(int i=0; i < self->nIds; ++i) {
        REQUIRE((u & self->tid[i]) == 0, file, lineno);
        u |= self->tid[i];
    }
#endif    
}

void IdSet_print(IdSet *self, FILE *fp) {
    if(self==NULL) {
        fputs("-----------------------\n", fp);
        return;
    }
    fprintf(fp, "pr=%lf nIds=%d:",
            self->p, self->nIds);
    for(int i=0; i < self->nIds; ++i)
        fprintf(fp, " 0%o", self->tid[i]);
    putc('\n', fp);
    fprintf(fp,"mig history:");
    MigOutcome_print(self->mig, fp);
}

/**
 * Allocate a new IdSet object with given values of nIds, tid, and
 * prob. Uses the struct hack.
 */
IdSet *IdSet_new(int nIds, const tipId_t *tid, double prob) {
    size_t size = sizeof(IdSet);
    if(nIds > 1)
        size += (nIds-1) * sizeof(tipId_t);
    
    IdSet *self = malloc(size);
    CHECKMEM(self);

    self->nIds = nIds;
    self->p = prob;
    self->mig = NULL;

    memcpy(self->tid, tid, nIds * sizeof(tipId_t));

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

/**
 * Join two IdSet objects, left and right, along with all the tipId_t
 * values in the "samples" array. On success, function returns a
 * pointer to the new IdSet object. If "left" and "right" represent
 * mutually exclusive events, they cannot be joined, and the function
 * returns NULL.
 */
IdSet *IdSet_join(IdSet *left, IdSet *right, int nsamples,
                  tipId_t *samples) {

    int mutually_exclusive;
    MigOutcome *mig = MigOutcome_join(left->mig, right->mig,
                                      &mutually_exclusive);
    if(mutually_exclusive)
        return NULL;

    int nIds = left->nIds + right->nIds + nsamples;
    
    // Copy all tipId_t values into tid. Allocating 1 extra
    // position to avoid problems with 0-length arrays.
    tipId_t tid[1+nIds];

    memcpy(tid, left->tid, left->nIds * sizeof(tipId_t));
    memcpy(tid + left->nIds, right->tid, right->nIds * sizeof(tipId_t));
    memcpy(tid + left->nIds + right->nIds, samples,
           nsamples * sizeof(tipId_t));

    IdSet *new = IdSet_new(nIds, tid, left->p * right->p);
    new->mig = mig;

    return new;
}

void IdSet_addMigEvent(IdSet *self, unsigned event, unsigned outcome,
                       double pr) {
    self->mig = MigOutcome_insert(self->mig, event, outcome, pr);
}

IdSet *IdSet_addSamples(IdSet *old, int nsamples, tipId_t *samples) {
    if(nsamples == 0)
        return old;
    
    int nIds = old->nIds + nsamples;
    tipId_t tid[nIds];

    memcpy(tid, old->tid, old->nIds * sizeof(tipId_t));
    memcpy(tid + old->nIds, samples, nsamples * sizeof(tipId_t));
    
    IdSet *new = IdSet_new(nIds, tid, old->p);
    new->mig = old->mig;
    old->mig = NULL;

    IdSet_free(old);
    return new;
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
