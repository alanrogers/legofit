#ifndef ARR_PARKEYVAL_H
#define ARR_PARKEYVAL_H

#include "typedefs.h"

void        ParKeyVal_free(ParKeyVal *node);
ParKeyVal  *ParKeyVal_add(ParKeyVal *node, const char *key, double *vptr);
int         ParKeyVal_get(ParKeyVal *node, double **valPtr, const char *key);

#endif
