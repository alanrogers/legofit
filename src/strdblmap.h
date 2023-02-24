#ifndef STRDBLMAP_H
#  define STRDBLMAP_H
#  include "typedefs.h"

StrDblMap *StrDblMap_new(int dim);
void       StrDblMap_free(StrDblMap * self);
double     StrDblMap_get(StrDblMap * self, const char * key, int *status);
int        StrDblMap_insert(StrDblMap * self, const char * key,
                                double value);
unsigned long StrDblMap_size(StrDblMap * self);
int        StrDblMap_hasKey(StrDblMap *self, const char *key);
int        StrDblMap_keys(StrDblMap *self, unsigned size, char * keys[size]);

#endif
