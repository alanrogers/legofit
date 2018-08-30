#ifndef LEGOFIT_STRPARMAP_H
#define LEGOFIT_STRPARMAP_H

#include "typedefs.h"

Param      *METHOD(search)(CLASS *root, KEYTYPE key);
CLASS      *METHOD(insert)(CLASS *root, Param *par);
void        METHOD(free)(CLASS *h);

#endif
