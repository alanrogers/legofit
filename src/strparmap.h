#ifndef LEGOFIT_STRPARMAP_H
#define LEGOFIT_STRPARMAP_H

#include "typedefs.h"

StrParMap  *StrParMap_insert(StrParMap *root, Param *par);
Param      *StrParMap_search(StrParMap *root, const char *key);
void        StrParMap_free(StrParMap *h);

#endif
