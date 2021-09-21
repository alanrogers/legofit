#ifndef ARR_LONGVEC_H
#define ARR_LONGVEC_H

#include "typedefs.h"

LongVec *LongVec_new(int size);
int      LongVec_size(const LongVec *self);
int      LongVec_resize(LongVec *self, int newsize);
void     LongVec_set(LongVec *self, int ndx, long value);
void     LongVec_plusEquals(LongVec *self, int ndx, long inc);
void     LongVec_minusEquals(LongVec *self, int ndx, long dec);
long     LongVec_get(const LongVec *self, int ndx);
void     LongVec_free(LongVec *self);

#endif
