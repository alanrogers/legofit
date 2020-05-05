#ifndef ARR_SAMPNDX_H
#  define ARR_SAMPNDX_H

#  include "typedefs.h"
#  include <stdlib.h>

struct SampNdx {
    // Array "node" contains an entry for each sample. That entry
    // is a pointer to the node into which the sample should
    // be placed. The sample gets a label of type tipIt_t. For sample
    // i, the label equals 2^i (i.e. 1<<i). There is another class,
    // called LblNdx, which maintains an array of labels. In that
    // array, the i'th label refers to the i'th sample in SampNdx.
    // I keep them separate, because LblNdx needs to be passed to
    // functions that have no need to know about pointers to PopNode
    // objects.
    unsigned    n;              // number of samples
    PopNode    *node[MAXSAMP];
};

void        SampNdx_init(SampNdx * self);
void        SampNdx_addSamples(SampNdx * self, unsigned nsamples,
							   PopNode * pnode);
void        SampNdx_populateTree(SampNdx * self);
unsigned    SampNdx_size(SampNdx * self);
int         SampNdx_equals(const SampNdx *lhs, const SampNdx *rhs);
void        SampNdx_sanityCheck(SampNdx *self, const char *file, int line);
int         SampNdx_ptrsLegal(SampNdx *self, PopNode *start, PopNode *end);
void        SampNdx_shiftPtrs(SampNdx *self, size_t dpop, int sign);

#endif
