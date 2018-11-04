#ifndef ARR_UINTQUEUE_H
#define ARR_UINTQUEUE_H

#include "typedefs.h"

UINTqueue *UINTqueue_push(UINTqueue *prev, unsigned val);
UINTqueue *UINTqueue_pop(UINTqueue *self, unsigned *value);
UINTqueue *UINTqueue_free(UINTqueue *self);
int UINTqueue_length(UINTqueue *self);

#endif
