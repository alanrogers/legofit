#ifndef STRSTRMAP_H
#  define STRSTRMAP_H
#  include "typedefs.h"

StrStrMap *StrStrMap_new(int dim);
void       StrStrMap_free(StrStrMap * self);
const char *StrStrMap_get(StrStrMap * self, const char * key, int *status);
int        StrStrMap_insert(StrStrMap * self, const char * key,
                                const char * value);
unsigned long StrStrMap_size(StrStrMap * self);
int        StrStrMap_hasKey(StrStrMap *self, const char *key);
int        StrStrMap_keys(StrStrMap *self, unsigned size, char * keys[size]);

#endif
