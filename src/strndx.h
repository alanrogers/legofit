#ifndef LEGO_STRNDX
#  define LEGO_STRNDX

#  include "typedefs.h"

StrNdx     *StrNdx_new(void);
void        StrNdx_free(StrNdx * self);
int         StrNdx_getNdx(StrNdx * self, char *key);
void        StrNdx_print(const StrNdx * self);
int         StrNdx_size(StrNdx * self);
#endif
