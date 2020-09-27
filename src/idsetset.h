#ifndef IDSETSET_H
#define IDSETSET_H

#include "typedefs.h"

int       IdSetSet_add(IdSetSet * self, IdSet *idset);
void      IdSetSet_empty_deep(IdSetSet * self);
void      IdSetSet_empty_shallow(IdSetSet * self);
void      IdSetSet_free_deep(IdSetSet * self);
void      IdSetSet_free_shallow(IdSetSet * self);
IdSetSet *IdSetSet_new(int dim);
IdSet    *IdSetSet_next(IdSetSet *self);
int       IdSetSet_reserve(IdSetSet *self, int m);
int       IdSetSet_rewind(IdSetSet *self);
void      IdSetSet_sanityCheck(IdSetSet *self, const char *file, int line);
int       IdSetSet_size(IdSetSet * self);
int       IdSetSet_toArray(IdSetSet *self, unsigned size, IdSet *v[size]);

#endif
