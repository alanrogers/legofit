#ifndef IDSET_TBL_H
#define IDSET_TBL_H

#include "typedefs.h"

IdSetTbl *IdSetTbl_new(int dim);
void      IdSetTbl_free(IdSetTbl * self);
int       IdSetTbl_add(IdSetTbl * self, IdSet *idset);
int       IdSetTbl_size(IdSetTbl * self);
int       IdSetTbl_toArray(IdSetTbl *self, unsigned size, IdSet *v[size]);
int       IdSetTbl_rewind(IdSetTbl *self);
IdSet    *IdSetTbl_next(IdSetTbl *self);

#endif
