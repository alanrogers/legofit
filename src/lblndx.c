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
const char *LblNdx_lbl(const LblNdx * self, unsigned i) {
    assert(i < self->n);
    return self->lbl[i];
}

unsigned LblNdx_size(const LblNdx * self) {
    return self->n;
}

void        LblNdx_sanityCheck(const LblNdx *self, const char *file, int line) {
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

int         LblNdx_equals(const LblNdx *lhs, const LblNdx *rhs) {
    if(lhs == rhs)
        return 1;
    if(lhs==NULL || rhs==NULL)
        return 0;
    if(lhs->n != rhs->n)
        return 0;

    int i;
    for(i=0; i < lhs->n; ++i)
        if(0 != strcmp(lhs->lbl[i], rhs->lbl[i]))
            return 0;
    return 1;
}

/// Reverse lookup. Return tipId_t value corresponding to
/// label. If label is not present in LblNdx, return 0.
tipId_t     LblNdx_getTipId(const LblNdx *self, const char *lbl) {
    unsigned i;
    tipId_t rval = 1;

    for(i=0; i < self->n; ++i) {
        if(0 == strcmp(lbl, self->lbl[i]))
           break;
    }
    if(i == self->n)
        return 0;

    rval <<= i;
    return rval;
}

void LblNdx_print(const LblNdx *self, FILE *fp) {
    unsigned i;
    for(i=0; i < self->n; ++i)
        fprintf(fp,"%3u %5s\n", i, self->lbl[i]);
}

/// Generate a label for site pattern tid. Label goes into
/// buff. Function returns a pointer to buff;
char       *patLbl(size_t n, char buff[n], tipId_t tid, const LblNdx * lblndx) {
    int         maxbits = 40;
    int         bit[maxbits];
    int         i, nbits;
    char        lbl[100];

    nbits = getBits(tid, maxbits, bit);
    buff[0] = '\0';
    for(i = 0; i < nbits; ++i) {
        snprintf(lbl, sizeof(lbl), "%s",
                 LblNdx_lbl(lblndx, (unsigned) bit[i]));
        if(strlen(buff) + strlen(lbl) >= n)
            eprintf("%s:%s:%d: buffer overflow\n", __FILE__, __func__,
                    __LINE__);
        strcat(buff, lbl);
        if(i + 1 < nbits && 1 + strlen(buff) < n)
            strcat(buff, ":");
    }
    return buff;
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

    for(i=0; i < LblNdx_size(&lndx); ++i) {
        assert((1u << i) == LblNdx_getTipId(&lndx, LblNdx_lbl(&lndx, i)));
    }

    if(verbose)
        LblNdx_print(&lndx, stdout);

    LblNdx lndx2 = {.n = 3 };
    LblNdx_init(&lndx2);
    LblNdx_addSamples(&lndx2, 1, "A");
    LblNdx_addSamples(&lndx2, 2, "B");

    assert(LblNdx_equals(&lndx, &lndx2));

	unitTstResult("LblNdx", "OK");

    return 0;
}
#endif
