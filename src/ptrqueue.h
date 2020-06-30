#ifndef PTRQUEUE_H
#define PTRQUEUE_H

#include "typedefs.h"

PtrQueue *PtrQueue_new(void);
void      PtrQueue_free(PtrQueue *self);
void      PtrQueue_push(PtrQueue *self, void *ptr);
void     *PtrQueue_pop(PtrQueue *self);
unsigned  PtrQueue_size(PtrQueue *self);

#endif
