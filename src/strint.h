#ifndef LEGO_STRINT
#  define LEGO_STRINT

#  include "typedefs.h"
#  include <stdio.h>

StrInt     *StrInt_new(void);
void        StrInt_free(StrInt * self);
void        StrInt_insert(StrInt *self, const char *key, int value);
int         StrInt_get(StrInt * self, const char *key);
void        StrInt_print(const StrInt * self, FILE *fp);
unsigned    StrInt_size(const StrInt *self);
#endif
