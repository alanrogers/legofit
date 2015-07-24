#ifndef HAVE_SAMPNDX
#  define HAVE_SAMPNDX

#  include "gptree.h"

struct SampNdx {
    unsigned    n;              // number of samples
    char        lbl[MAXSAMP][POPNAMESIZE];
    PopNode    *node[MAXSAMP];
};

void        SampNdx_init(SampNdx * self);
void        SampNdx_addSamples(SampNdx * self, unsigned nsamples,
                               PopNode * pnode, const char *lbl);
void        SampNdx_populateTree(SampNdx * self);
const char *SampNdx_lbl(SampNdx * self, unsigned i);
unsigned    SampNdx_size(SampNdx * self);

#endif
