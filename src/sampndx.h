#ifndef ARR_SAMPNDX_H
#  define ARR_SAMPNDX_H

#  include "typedefs.h"
#  include <stdlib.h>

struct SampNdx {
    // Array "ptr" contains an entry for each sample. That entry is a
    // pointer to "home", the node or segment that contains the
    // sample.  For sample i, the label equals 2^i (i.e. 1<<i).
    unsigned  n;              // number of samples
    void     *ptr[MAXSAMP];   // array of pointers
};

void        SampNdx_init(SampNdx * self);
void        SampNdx_addSamples(SampNdx * self, unsigned nsamples, void * ptr);
void       *SampNdx_get(SampNdx * self, unsigned i);
unsigned    SampNdx_size(SampNdx * self);
int         SampNdx_equals(const SampNdx *lhs, const SampNdx *rhs);
void        SampNdx_sanityCheck(SampNdx *self, const char *file, int line);
int         SampNdx_ptrsLegal(SampNdx * self, void * start, void * end);
<<<<<<< HEAD
void        SampNdx_shiftPtrs(SampNdx *self, size_t dpop, int sign);
=======
void        SampNdx_remapPtrs(SampNdx *self, PtrPtrMap *ppm);
>>>>>>> matcoal

#endif
