#ifndef PTRPTRMAP_H
#  define PTRPTRMAP_H
#  include "typedefs.h"

PtrPtrMap *PtrPtrMap_new(void);
void       PtrPtrMap_free(PtrPtrMap * self);
void      *PtrPtrMap_get(PtrPtrMap * self, const void *key);
int        PtrPtrMap_insert(PtrPtrMap * self, const void *key, void * value);
unsigned long PtrPtrMap_size(PtrPtrMap * self);
void       PtrPtrMap_print(PtrPtrMap * self);

#endif
