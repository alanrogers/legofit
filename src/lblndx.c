/**
 * @file lblndx.c
 * @brief An index of sample labels.
 */
#include "lblndx.h"
#include "misc.h"
#include "parkeyval.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

/// Set everything to zero.
void LblNdx_init(LblNdx * self) {
    memset(self, 0, sizeof(*self));
    LblNdx_sanityCheck(self, __FILE__, __LINE__);
}

/// Add samples for a single population. Should be called once for
/// each sampled population.
void LblNdx_addSamples(LblNdx * self, unsigned nsamples, const char *lbl) {
    unsigned    i;
    if(self->n + nsamples >= MAXSAMP)
        eprintf("%s:%s:%d: too many samples\n", __FILE__, __func__, __LINE__);
    for(i = 0; i < nsamples; ++i) {
        if(nsamples == 1)
            snprintf(self->lbl[self->n], POPNAMESIZE, "%s", lbl);
        else
            snprintf(self->lbl[self->n], POPNAMESIZE, "%s.%u", lbl, i);
        self->n += 1;
    }
    LblNdx_sanityCheck(self, __FILE__, __LINE__);
}

/// Return the label associated with index i.
const char *LblNdx_lbl(LblNdx * self, unsigned i) {
    assert(i < self->n);
    return self->lbl[i];
}

unsigned LblNdx_size(LblNdx * self) {
    return self->n;
}

void        LblNdx_sanityCheck(LblNdx *self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self, file, line);
    REQUIRE(self->n < MAXSAMP, file, line);

    int i;
    for(i=0; i < self->n; ++i) {
        REQUIRE(NULL != self->lbl[i], file, line);
        REQUIRE(legalName(self->lbl[i]), file, line);
    }
#endif
}

#ifdef TEST

#  include <string.h>
#  include <assert.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

int main(int argc, char **argv) {

    int         verbose = 0;
    unsigned    i;

    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xlblndx [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    LblNdx     lndx = {.n = 3 };
    assert(lndx.n == 3);
    assert(LblNdx_size(&lndx) == 3);

    LblNdx_init(&lndx);
    assert(LblNdx_size(&lndx) == 0);

    LblNdx_addSamples(&lndx, 1, "A");
    LblNdx_addSamples(&lndx, 2, "B");

    assert(3 == LblNdx_size(&lndx));
    assert(0 == strcmp("A", LblNdx_lbl(&lndx, 0)));
    assert(0 == strcmp("B.0", LblNdx_lbl(&lndx, 1)));
    assert(0 == strcmp("B.1", LblNdx_lbl(&lndx, 2)));

    for(i = 0; verbose && i < LblNdx_size(&lndx); ++i) {
        printf("%d %s\n", i, LblNdx_lbl(&lndx, i));
    }

	unitTstResult("LblNdx", "OK");

    return 0;
}
#endif
