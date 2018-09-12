#ifndef LEGOFIT_PTRSET_H
#define LEGOFIT_PTRSET_H

#include "typedefs.h"

PtrSet *PtrSet_new(void);
void    PtrSet_free(PtrSet *self);
void    PtrSet_insert(PtrSet *self, double *ptr);
int     PtrSet_exists(PtrSet *self, double *ptr);

#endif
