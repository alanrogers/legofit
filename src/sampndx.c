/**
 * @file sampndx.c
 * @brief An index of samples.
 */
#include "sampndx.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

/// Set everything to zero.
void SampNdx_init(SampNdx * self) {
    memset(self, 0, sizeof(*self));
}

/// Add samples for a single population. Should be called once for
/// each sampled population.
void SampNdx_addSamples(SampNdx * self, unsigned nsamples, PopNode * pnode,
                        const char *lbl) {
    unsigned    i;
    if(self->n + nsamples >= MAXSAMP)
        eprintf("%s:%s:%d: too many samples\n", __FILE__, __func__, __LINE__);
    for(i = 0; i < nsamples; ++i) {
        if(nsamples == 1)
            snprintf(self->lbl[self->n], POPNAMESIZE, "%s", lbl);
        else
            snprintf(self->lbl[self->n], POPNAMESIZE, "%s.%u", lbl, i);
        self->node[self->n] = pnode;
        self->n += 1;
    }
}

/// Put samples into the gene tree. Should be done at the start of
/// each simulation.
void SampNdx_populateTree(SampNdx * self) {
    unsigned    i;
    for(i = 0; i < self->n; ++i)
        PopNode_newGene(self->node[i], i);
}

/// Return the label associated with index i.
const char *SampNdx_lbl(SampNdx * self, unsigned i) {
    assert(i < self->n);
    return self->lbl[i];
}

unsigned SampNdx_size(SampNdx * self) {
    return self->n;
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

    PopNode    *pnode = PopNode_new(100.0, 20.0);
    SampNdx_addSamples(&sndx, 1, pnode, "A");
    SampNdx_addSamples(&sndx, 2, pnode, "B");

    assert(3 == SampNdx_size(&sndx));
    assert(0 == strcmp("A", SampNdx_lbl(&sndx, 0)));
    assert(0 == strcmp("B.0", SampNdx_lbl(&sndx, 1)));
    assert(0 == strcmp("B.1", SampNdx_lbl(&sndx, 2)));

    for(i = 0; verbose && i < SampNdx_size(&sndx); ++i) {
        printf("%d %s\n", i, SampNdx_lbl(&sndx, i));
    }

    SampNdx_populateTree(&sndx);
    assert(3 == PopNode_nsamples(pnode));

    return 0;
}
#endif
