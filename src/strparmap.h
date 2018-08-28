#ifndef LEGOFIT_STRPARMAP_H
#define LEGOFIT_STRPARMAP_H

#include "typedefs.h"

const char *StrParMap_key(StrParMap *h);
Param      *StrParMap_search(StrParMap *root, const char *key);
StrParMap  *StrParMap_insert(StrParMap *root, Param *par);
void        StrParMap_free(StrParMap *h);

#endif
