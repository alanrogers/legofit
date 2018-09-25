#ifndef LEGOFIT_ADDRPARMAP_H
#define LEGOFIT_ADDRPARMAP_H
#include "typedefs.h"

AddrParMap *AddrParMap_insert(AddrParMap *root, Param *par);
Param      *AddrParMap_search(AddrParMap *root, const double *valptr);
void        AddrParMap_free(AddrParMap *h);

#endif
