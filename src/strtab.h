#ifndef LEGO_STRTAB
#  define LEGO_STRTAB

#  include "typedefs.h"

StrTab     *StrTab_new(void);
void        StrTab_free(StrTab * self);
int         StrTab_get(StrTab * self, char *key);
void        StrTab_add(StrTab * self, char *key, int value);
unsigned    StrTab_size(StrTab * self);
void        StrTab_print(const StrTab * self);
#endif
