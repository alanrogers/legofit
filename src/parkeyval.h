#ifndef ARR_PARKEYVAL_H
#define ARR_PARKEYVAL_H

#include "typedefs.h"
#include <stdio.h>

void        ParKeyVal_free(ParKeyVal *node);
ParKeyVal  *ParKeyVal_add(ParKeyVal *node, const char *key, double *vptr);
double     *ParKeyVal_get(ParKeyVal *node, const char *key);
void        ParKeyVal_print(ParKeyVal *self, FILE *fp);

#endif
