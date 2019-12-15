#ifndef ARR_SETPART_H
#define ARR_SETPART_H

#include "stirling2.h"

typedef struct SetPart SetPart;

SetPart      *SetPart_new(unsigned nelements, unsigned nsubdivisions,
                          Stirling2 *s);
void          SetPart_free(SetPart *self);
unsigned      SetPart_sizeOfSet(SetPart *self);
unsigned      SetPart_nSubsets(SetPart *self);
long unsigned SetPart_nPartitions(SetPart *self);
void          SetPart_getPartition(const SetPart *self, long unsigned i,
                                   unsigned n, unsigned p[n]);
#endif
