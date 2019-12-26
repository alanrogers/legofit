#ifndef ARR_SEGMENT_H
#define ARR_SEGMENT_H

typedef struct Segment Segment;

int Segment_coalesce(Segment *self, int maxsamp,
                     MatCoal *mc[maxsamp-1],
                     double v);

#endif
