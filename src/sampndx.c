/**
 * @file sampndx.c
 * @author Alan R. Rogers
 * @brief An index of haploid genetic samples
 *
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "sampndx.h"
#include "misc.h"
#include <string.h>
#include <stdio.h>
/// Set everything to zero.
void SampNdx_init(SampNdx * self) {
    memset(self, 0, sizeof(*self));
}

/// Add samples for a single population. Should be called once for
/// each sampled population. This justs sets a pointer and increments
/// the count of pointers. It doesn't allocate anything.
void SampNdx_addSamples(SampNdx * self, unsigned nsamples, void * ptr) {
    unsigned    i;
    if(self->n + nsamples >= MAXSAMP)
        eprintf("%s:%s:%d: too many samples\n", __FILE__, __func__, __LINE__);
    for(i = 0; i < nsamples; ++i) {
        self->ptr[self->n] = ptr;
        self->n += 1;
    }
}

/// Return pointer to home of i'th sample.
void *SampNdx_get(SampNdx * self, unsigned i) {
    assert(i < self->n);
    return self->ptr[i];
}

unsigned SampNdx_size(SampNdx * self) {
    return self->n;
}

/// This equality check doesn't do much, because the pointers in
/// different SampNdx objects don't have to be (in fact shouldn't be)
/// equal.
int SampNdx_equals(const SampNdx * lhs, const SampNdx * rhs) {
    if(lhs == NULL && rhs == NULL)
        return 1;
    if(lhs == NULL || rhs == NULL)
        return 0;
    if(lhs->n != rhs->n)
        return 0;
    return 1;
}

/// Check sanity of a SampNdx.
void SampNdx_sanityCheck(SampNdx * self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self != NULL, file, line);
    REQUIRE(self->n < MAXSAMP, file, line);
    for(int i = 0; i < self->n; ++i)
        REQUIRE(NULL != self->ptr[i], file, line);
#endif
}

/// Return 1 if all pointers in SampNdx are in [start,end); return 0
/// otherwise.
int SampNdx_ptrsLegal(SampNdx * self, void * start, void * end) {
    assert(self);
    for(int i = 0; i < self->n; ++i) {
        if(self->ptr[i] < start || self->ptr[i] >= end)
            return 0;
    }
    return 1;
}

/// Shift all pointers within SampNdx by an offset of magnitude
/// dpop. If sign > 0, the shift is positive; otherwise it is
/// negative.
void SampNdx_shiftPtrs(SampNdx * self, size_t dpop, int sign) {
    for(int i = 0; i < self->n; ++i)
        SHIFT_PTR(self->ptr[i], dpop, sign);
}

#ifdef TEST

#  include <string.h>
#  include <assert.h>
#  include <time.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

int main(int argc, char **argv) {

    int         verbose = 0;

    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xsampndx [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    SampNdx     sndx = {.n = 3 };
    assert(sndx.n == 3);
    assert(SampNdx_size(&sndx) == 3);

    SampNdx_init(&sndx);
    assert(SampNdx_size(&sndx) == 0);

    int         nseg = 3;
    int         v[nseg];

    // Each sample has a different pointer, although this isn't
    // necessary.
    SampNdx_addSamples(&sndx, 1, &v[0]);
    SampNdx_addSamples(&sndx, 1, &v[1]);
    SampNdx_addSamples(&sndx, 1, &v[2]);

    assert(3 == SampNdx_size(&sndx));

    // Make sure all the pointers are within [v, v+nseg)
    assert(SampNdx_ptrsLegal(&sndx, v, v + nseg));

    assert(SampNdx_equals(&sndx, &sndx));

    for(int i=0; i<3; ++i)
        assert(&v[i] == SampNdx_get(&sndx, i));

    SampNdx_sanityCheck(&sndx, __FILE__, __LINE__);

    SampNdx sndx2 = sndx;

    SampNdx_addSamples(&sndx2, 1, &v[0]);
    assert(4 == SampNdx_size(&sndx2));
    assert(!SampNdx_equals(&sndx, &sndx2));


    unitTstResult("SampNdx", "OK");

    return 0;
}
#endif
    
