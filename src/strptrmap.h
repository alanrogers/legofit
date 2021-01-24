#ifndef STRPTRMAP_H
#  define STRPTRMAP_H
#  include "typedefs.h"

StrPtrMap *StrPtrMap_new(void);
void       StrPtrMap_free(StrPtrMap * self);
void      *StrPtrMap_get(StrPtrMap * self, const char *key);
int        StrPtrMap_insert(StrPtrMap * self, const char *key,
                            void * node);
unsigned long StrPtrMap_size(StrPtrMap * self);
void       StrPtrMap_print(StrPtrMap * self);
void       StrPtrMap_ptrArray(StrPtrMap *self, long unsigned n, void *v[n]);

#endif
