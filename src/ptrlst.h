#ifndef ARR_PTRLST_H
#define ARR_PTRLST_H

#include "typedefs.h"
#include <assert.h>

typedef struct PtrLstEl PtrLstEl;

PtrLst *PtrLst_new(void);
void    PtrLst_free(PtrLst *self);
int     PtrLst_push(PtrLst *self, void *ptr);
void   *PtrLst_pop(PtrLst *self);
long unsigned PtrLst_length(PtrLst *self);
void    PtrLst_rewind(PtrLst *self);
void   *PtrLst_next(PtrLst *self);

#endif
