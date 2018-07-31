#ifndef LEGO_STRDBLMAP
#  define LEGO_STRDBLMAP

#  include "typedefs.h"
#  include <stdio.h>

StrDblMap  *StrDblMap_new(void);
void        StrDblMap_free(StrDblMap * self);
void        StrDblMap_insert(StrDblMap *self, const char *key, double value);
double      StrDblMap_get(StrDblMap * self, const char *key);
int         StrDblMap_exists(StrDblMap * self, const char *key);
void        StrDblMap_print(const StrDblMap * self, FILE *fp);
unsigned    StrDblMap_size(const StrDblMap *self);
#endif
