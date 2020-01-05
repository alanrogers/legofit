#ifndef ARR_SEGMENT_H
#define ARR_SEGMENT_H

#include "typedefs.h"

typedef struct Segment Segment;

int Segment_coalesce(Segment *self, int maxsamp, BranchTab *branchtab,
                     double v);

#endif
