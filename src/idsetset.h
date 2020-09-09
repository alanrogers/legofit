#ifndef IDSETSET_H
#define IDSETSET_H

#include "typedefs.h"

IdSetSet *IdSetSet_new(int dim);
void      IdSetSet_free(IdSetSet * self);
int       IdSetSet_add(IdSetSet * self, IdSet *idset);
int       IdSetSet_size(IdSetSet * self);
int       IdSetSet_toArray(IdSetSet *self, unsigned size, IdSet *v[size]);
int       IdSetSet_rewind(IdSetSet *self);
IdSet    *IdSetSet_next(IdSetSet *self);

#endif
