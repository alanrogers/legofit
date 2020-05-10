#ifndef ARR_SAMPNDX_H
#  define ARR_SAMPNDX_H

#  include "typedefs.h"
#  include <stdlib.h>

struct SampNdx {
    // Array "node" contains an entry for each sample. That entry is a
    // pointer to "home", the node or segment that contains the
    // sample.  This is a void pointer, so that SampNdx doesn't have
    // to know whether it's storing a pointer to a PopNode or to a
    // Segment.  For sample i, the label equals 2^i (i.e. 1<<i). There
    // is another class, called LblNdx, which maintains an array of
    // labels. In that array, the i'th label refers to the i'th sample
    // in SampNdx.  I keep them separate, because LblNdx needs to be
    // passed to functions that have no need to know about pointers to
    // PopNode objects.
    unsigned  n;              // number of samples
    void     *ptr[MAXSAMP];  // array of pointers
};

void        SampNdx_init(SampNdx * self);
void        SampNdx_addSamples(SampNdx * self, unsigned nsamples, void * ptr);
void       *SampNdx_get(SampNdx * self, unsigned i);
unsigned    SampNdx_size(SampNdx * self);
int         SampNdx_equals(const SampNdx *lhs, const SampNdx *rhs);
void        SampNdx_sanityCheck(SampNdx *self, const char *file, int line);
int         SampNdx_ptrsLegal(SampNdx * self, void * start, void * end);
void        SampNdx_shiftPtrs(SampNdx *self, size_t dpop, int sign);

#endif
