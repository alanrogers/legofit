#ifndef LEGOFIT_PTRSET_H
#define LEGOFIT_PTRSET_H

#include "typedefs.h"
#include <stdio.h>

PtrSet *PtrSet_new(void);
void    PtrSet_free(PtrSet *self);
void    PtrSet_insert(PtrSet *self, const double *ptr);
int     PtrSet_exists(PtrSet *self, const double *ptr);
void    PtrSet_print(PtrSet *self, FILE *fp);

#endif
