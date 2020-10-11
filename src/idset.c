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

static void merge(int nz, tipId_t *z, int nx, tipId_t *x,
                  int ny, tipId_t *y);

void IdSet_sanityCheck(IdSet *self, const char *file, int lineno) {
#ifndef NDEBUG
    if(self == NULL)
       return;
    REQUIRE(self->nIds >= 0, file, lineno);
    REQUIRE(self->p >= 0.0, file, lineno);
    REQUIRE(self->p <= 1.0, file, lineno);
    for(int i=0; i < self->nIds; ++i)
        REQUIRE(self->tid[i] > 0, file, lineno);

    // Empty IdSet objects arise only by migration.
    if(self->nIds == 0)
        REQUIRE(NULL !=  self->evlst, file, lineno);

    // tipId_t values should be sorted in increasing order
    for(int i=1; i < self->nIds; ++i)
        REQUIRE(self->tid[i-1] < self->tid[i], file, lineno);

    // tipId_t values should not share bits
    REQUIRE(no_shared_bits(self->nIds, self->tid), file, lineno);
#endif    
}

void IdSet_print(IdSet *self, FILE *fp) {
#if 0    
    fprintf(fp, "pr=%Lf nIds=%d:",
            self->p, self->nIds);
#endif
    fputs("[", fp);
    for(int i=0; i < self->nIds; ++i) {
        if(i>0)
            putc(' ', fp);
        fprintf(fp, "%o", self->tid[i]);
    }
    fprintf(fp,":%Lg", self->p);
    putc(':', fp);
    Event_print(self->evlst, fp);
    fputs("]", fp);
}

/**
 * Allocate a new IdSet object with given values of nIds, tid, and
 * prob. Uses the struct hack.
 */
IdSet *IdSet_new(int nIds, const tipId_t *tid, long double prob) {
    size_t size = sizeof(IdSet);
    if(nIds > 1)
        size += (nIds-1) * sizeof(tipId_t);
    
    IdSet *self = malloc(size);
    CHECKMEM(self);

    self->nIds = nIds;
    self->p = prob;
    self->evlst = NULL;

    memcpy(self->tid, tid, nIds * sizeof(tipId_t));

    return self;
}

// Duplicate
IdSet *IdSet_dup(const IdSet *old) {
    if(old == NULL)
        return NULL;
    IdSet *new = IdSet_new(old->nIds, old->tid, old->p);
    new->evlst = Event_dup(old->evlst);
    return new;
}

// Add Event list from old into self.
void IdSet_copyEvent(IdSet *self, const IdSet *old) {
    assert(self->evlst == NULL);
    self->evlst = Event_dup(old->evlst);
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
    Event *evlst = Event_join(left->evlst, right->evlst,
                              &mutually_exclusive);
    if(mutually_exclusive)
        return NULL;

    int left_plus_right = left->nIds + right->nIds;
    int nIds = left_plus_right + nsamples;
    
    // Copy all tipId_t values into tid. Allocating 1 extra
    // position to avoid problems with 0-length arrays.
    tipId_t tid[1+nIds], buff[1+left_plus_right];

    // Copy all ids into tid while maintaining sort.
    merge(left_plus_right, buff, left->nIds, left->tid,
          right->nIds, right->tid);
    merge(nIds, tid, left_plus_right, buff, nsamples, samples);

    IdSet *new = IdSet_new(nIds, tid, left->p * right->p);
    new->evlst = evlst;

#ifndef NDEBUG    
    if(!no_shared_bits(nIds, tid)) {
        fprintf(stderr,"%s:%s:%d: shared bits\n",
                __FILE__,__func__,__LINE__);

        putc('\n', stderr);
        fprintf(stderr,"new:");
        IdSet_print(new, stderr);

        putc('\n', stderr);
        fprintf(stderr,"left:");
        IdSet_print(left, stderr);

        putc('\n', stderr);
        fprintf(stderr,"right:");
        IdSet_print(right, stderr);
        putc('\n', stderr);
    }
#endif    

    return new;
}

void IdSet_addEvent(IdSet *self, unsigned event, unsigned outcome,
                       long double pr) {
    self->event = Event_insert(self->event, event, outcome, pr);
}

/// Return pointer to new IdSet and free old one.
IdSet *IdSet_addSamples(IdSet *old, int nsamples, tipId_t *samples) {
    if(nsamples == 0)
        return old;
    
    int nIds = old->nIds + nsamples;
    tipId_t tid[nIds];

    merge(nIds, tid, old->nIds, old->tid, nsamples, samples);

    IdSet *new = IdSet_new(nIds, tid, old->p);
    new->evlst = old->evlst;
    old->evlst = NULL;

    IdSet_free(old);
    return new;
}

/// Allocate a new IdSet with a single tipId_t value and probability 1.
IdSet *IdSet_newTip(tipId_t tid) {
    IdSet *self = malloc(sizeof(IdSet));
    CHECKMEM(self);

    self->nIds = 1;
    self->p = 1.0;
    self->evlst = NULL;
    self->tid[0] = tid;
    return self;
}

int IdSet_cmp(const IdSet *left, const IdSet *right) {
    int diff = left->nIds - right->nIds;
    if(diff)
        return diff;
    for(int i=0; i < left->nIds; ++i) {
        if(left->tid[i] > right->tid[i])
            return 1;
        else if(left->tid[i] < right->tid[i])
            return -1;
    }
    return Event_cmp(left->evlst, right->evlst);
}

uint32_t IdSet_hash(const IdSet *self)
    __attribute__((no_sanitize("integer"))) {
    uint32_t hash = 17;
    for(int i=0; i < self->nIds; ++i)
        hash = hash * 37 + uint32Hash(self->tid[i]);
    uint32_t h = Event_hash(self->evlst);
    if(h)
        hash = hash * 37 + h;
    return hash;
}

// Copy x and y into z, while maintaining sort.
static void merge(int nz, tipId_t *z, int nx, tipId_t *x,
                  int ny, tipId_t *y) {
    assert(nz == nx + ny);

    int ix=0, iy=0, iz=0;

    while(ix<nx && iy<ny) {
        if(x[ix] < y[iy])
            z[iz++] = x[ix++];
        else
            z[iz++] = y[iy++];
    }
    while(ix < nx)
        z[iz++] = x[ix++];
    while(iy < ny)
        z[iz++] = y[iy++];
}
