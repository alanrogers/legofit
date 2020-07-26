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

void IdSet_sanityCheck(IdSet *self, const char *file, int lineno) {
#ifndef NDEBUG
    if(self == NULL)
       return;
    REQUIRE(self->nIds > 0, file, lineno);
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

/**
 * Join two IdSet objects, left and right, along with all the tipId_t
 * values in the "samples" array. On success, function returns a
 * pointer to the new IdSet object. If "left" and "right" represent
 * mutually exclusive events, they cannot be joined, and the function
 * returns NULL.
 */
IdSet *IdSet_join(IdSet *left, IdSet *right, int nsamples,
                  tipId_t *samples) {

    MigOutcome *mig = MigOutcome_join(left->mig, right->mig);
    if(mig == NULL)
        return NULL; // left and right are mutually exclusive

    // Copy all tipId_t values into tid, excluding zeroes, which
    // represent the empty set.
    int nIds = left->nIds + right->nIds + nsamples;
    tipId_t tid[nIds];
    for(int i=nIds=0; i < left->nIds; ++i) {
        if(left->tid[i])
            tid[nIds++] = left->tid[i];
    }
    for(int i=0; i < right->nIds; ++i) {
        if(right->tid[i])
            tid[nIds++] = left->tid[i];
    }
    for(int i=0; i<nsamples; ++i)
        tid[nIds++] = samples[i];

    // In case left, right, and samples are all empty.
    if(nIds == 0) {
        nIds += 1;
        tid[0] = 0; // the empty set
    }

    IdSet *new = IdSet_new(nIds, tid, left->p * right->p);
    new->mig = mig;

    return new;
}

void IdSet_addMigEvent(IdSet *self, unsigned event, unsigned outcome,
                       double pr) {
    self->mig = MigOutcome_insert(self->mig, event, outcome, pr);
}

