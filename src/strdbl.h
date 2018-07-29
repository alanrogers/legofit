#ifndef LEGO_STRDBL
#  define LEGO_STRDBL

#  include "typedefs.h"
#  include <stdio.h>

StrDbl     *StrDbl_new(void);
void        StrDbl_free(StrDbl * self);
void        StrDbl_insert(StrDbl *self, const char *key, double value);
double      StrDbl_get(StrDbl * self, const char *key);
int         StrDbl_exists(StrDbl * self, const char *key);
void        StrDbl_print(const StrDbl * self, FILE *fp);
unsigned    StrDbl_size(const StrDbl *self);
#endif
