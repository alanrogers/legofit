#ifndef HASHTAB_INCLUDED
#  define HASHTAB_INCLUDED

#include "typedefs.h"

HashTab    *HashTab_new(void);
void *      El_get(El * self);
void        El_set(El * self, void * value);
void        HashTab_free(HashTab * self);
El         *HashTab_get(HashTab * self, const char *key);
unsigned long HashTab_size(HashTab * self);
void        HashTab_print(HashTab *self);
#endif
