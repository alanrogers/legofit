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

// A set of tipId_t values.
struct IdSet {
    int nIds; // number of Ids in this set

    EventLst *evlst; // linked list of events, outcomes and probs.

    // Array of length nIds. tid[i] is the id of the i'th haploid
    // individual in this set. Allocated using the struct hack.
    tipId_t tid[1];
};

static void merge(int nz, tipId_t *z, int nx, tipId_t *x,
                  int ny, tipId_t *y);
static void merge3(int nw, tipId_t *w, int nz, tipId_t *z,
                   int nx, tipId_t *x, int ny, tipId_t *y);

void IdSet_sanityCheck(IdSet *self, const char *file, int lineno) {
#ifndef NDEBUG
    if(self == NULL)
       return;
    REQUIRE(self->nIds >= 0, file, lineno);
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
    fprintf(fp, " nIds=%d:", self->nIds);
#endif
    fputs("[", fp);
    for(int i=0; i < self->nIds; ++i) {
        if(i>0)
            putc(' ', fp);
        fprintf(fp, "%u", self->tid[i]);
    }
    putc(':', fp);
    EventLst_print(self->evlst, fp);
    fputs("]", fp);
}

/**
 * Allocate a new IdSet object with given values of nIds, tid, and
 * prob. Uses the struct hack.
 */
IdSet *IdSet_new(int nIds, const tipId_t *tid, EventLst *evlst) {
    size_t size = sizeof(IdSet);
    if(nIds > 1)
        size += (nIds-1) * sizeof(tipId_t);
    
    IdSet *self = malloc(size);
    CHECKMEM(self);

    self->nIds = nIds;
    self->evlst = evlst;

    memcpy(self->tid, tid, nIds * sizeof(tipId_t));

    return self;
}

// Duplicate
IdSet *IdSet_dup(const IdSet *old) {
    if(old == NULL)
        return NULL;
    IdSet *new = IdSet_new(old->nIds, old->tid,
                           EventLst_dup(old->evlst));
    return new;
}

// Add Event list from old into self.
void IdSet_copyEventLst(IdSet *self, const IdSet *old) {
    assert(self->evlst == NULL);
    self->evlst = EventLst_dup(old->evlst);
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
    EventLst *evlst = EventLst_join(left->evlst, right->evlst,
                                    &mutually_exclusive);
    if(mutually_exclusive)
        return NULL;

    int nIds = left->nIds + right->nIds + nsamples;
    
    // Copy all tipId_t values into tid. Allocating 1 extra
    // position to avoid problems with 0-length arrays.
    tipId_t tid[1+nIds];

    // Copy all ids into tid while maintaining sort.
    if(nsamples) {
        merge3(nIds, tid,
               left->nIds, left->tid,
               right->nIds, right->tid,
               nsamples, samples);
    }else{
        merge(nIds, tid,
              left->nIds, left->tid,
              right->nIds, right->tid);
    }

    IdSet *new = IdSet_new(nIds, tid, evlst);

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
    self->evlst = EventLst_insert(self->evlst, event, outcome, pr);
}

/// Return pointer to new IdSet and free old one.
IdSet *IdSet_addSamples(IdSet *old, int nsamples, tipId_t *samples) {
    if(nsamples == 0)
        return old;
    
    int nIds = old->nIds + nsamples;
    tipId_t tid[nIds];

    merge(nIds, tid, old->nIds, old->tid, nsamples, samples);

    IdSet *new = IdSet_new(nIds, tid, old->evlst);
    old->evlst = NULL;

    IdSet_free(old);
    return new;
}

/// Allocate a new IdSet with a single tipId_t value and probability 1.
IdSet *IdSet_newTip(tipId_t tid) {
    IdSet *self = malloc(sizeof(IdSet));
    CHECKMEM(self);

    self->nIds = 1;
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
    return EventLst_cmp(left->evlst, right->evlst);
}

uint32_t IdSet_hash(const IdSet *self)
    __attribute__((no_sanitize("integer"))) {
    uint32_t hash = 17;
    for(int i=0; i < self->nIds; ++i)
        hash = hash * 37 + uint32Hash(self->tid[i]);
    uint32_t h = EventLst_hash(self->evlst);
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

// Copy x, y, and z into w while maintaining sort.
static void merge3(int nw, tipId_t *w, int nz, tipId_t *z,
                   int nx, tipId_t *x, int ny, tipId_t *y) {
    assert(nz == nx + ny + nz);

    int iw=0, ix=0, iy=0, iz=0;

    while(ix<nx && iy<ny && iz<nz) {
        if(x[ix] < y[iy]) {
            if(x[ix] < z[iz])
                w[iw++] = x[ix++];
            else
                w[iw++] = z[iz++];
        }else{
            // y <= x
            if(y[iy] < z[iz])
                w[iw++] = y[iy++];
            else
                w[iw++] = z[iz++];
        }
    }

    if(iz==nz) {
        merge(nw-iw, w+iw, nx-ix, x+ix, ny-iy, y+iy);
    }else if(iy==ny) {
        merge(nw-iw, w+iw, nx-ix, x+ix, nz-iz, z+iz);
    }else{
        // ix==nx
        merge(nw-iw, w+iw, ny-iy, y+iy, nz-iz, z+iz);
    }
}

void IdSet_free(IdSet *self) {
    EventLst_free(self->evlst);
    free(self);
}

long double IdSet_prob(IdSet *self) {
    return EventLst_prob(self->evlst);
}

int IdSet_nIds(IdSet *self) {
    return self->nIds;
}

/// Return a tipId_t value representing the union of the values
/// within this IdSet object whose indices are given in ndx, an array
/// of length n.
tipId_t IdSet_union(IdSet *self, int n, int *ndx) {
    tipId_t u = 0;
    for(int j=0; j < n; ++j) {
        assert(ndx[j] < self->nIds);
        assert(ndx[j] >= 0);
        assert(self->tid[ndx[j]]);
        u |= self->tid[ndx[j]];
    }
    return u;
}

/// On entry, k is the number of parts into which this IdSet will be
/// subdivided, n must equal self->nIds, and a[j] is the index of the
/// part into which the j'th tipId_t value will be assigned.
///
/// On return, partition[i] is the union of the tipId_t values in self
/// which belong to the i'th part.
void IdSet_partition(IdSet *self, unsigned k, tipId_t part[k],
                    unsigned n, unsigned a[n]) {
    memset(part, 0, k*sizeof(tipId_t));
    assert(n == self->nIds);
    
    // Loop over tid values, "or"ing each one into the relevant
    // entry of "part", as specified in array "a". 
    for(int j=0; j<n; ++j) {
        assert(a[j] < k);
        part[a[j]] |= self->tid[j];
    }
}

/// Return a duplicate of the EventLst of this IdSet.
EventLst *IdSet_dupEventLst(const IdSet *self) {
    return EventLst_dup(self->evlst);
}

/// Return a pointer to an IdSet representing a subset of "self".
/// The tipId_t values of the subset are defined by "ndx", an array
/// of n non-negative integers, which index values in the original
/// IdSet.
IdSet *IdSet_subset(const IdSet *self, int n, int *ndx) {
    tipId_t tid[n+1];
    for(int i=0; i<n; ++i) {
        assert(ndx[i] >= 0);
        assert(ndx[i] < self->nIds);
        tid[i] = self->tid[ndx[i]];
    }

    EventLst *evlst = EventLst_dup(self->evlst);
    return IdSet_new(n, tid, evlst);
}

/// Return i'th tipId_t value.
int    IdSet_get(const IdSet *self, int i) {
    assert(i>=0);
    assert(i < self->nIds);
    return(self->tid[i]);
}

