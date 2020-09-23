#ifndef IDSETSET_H
#define IDSETSET_H

#include "typedefs.h"

int       IdSetSet_add(IdSetSet * self, IdSet *idset);
void      IdSetSet_empty(IdSetSet * self);
void      IdSetSet_free(IdSetSet * self);
IdSetSet *IdSetSet_new(int dim);
IdSet    *IdSetSet_next(IdSetSet *self);
int       IdSetSet_rewind(IdSetSet *self);
int       IdSetSet_size(IdSetSet * self);
int       IdSetSet_toArray(IdSetSet *self, unsigned size, IdSet *v[size]);

#endif
