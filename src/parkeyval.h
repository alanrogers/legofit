#ifndef ARR_PARKEYVAL_H
#define ARR_PARKEYVAL_H

#include "typedefs.h"
#include <stdio.h>

void        ParKeyVal_free(ParKeyVal *node);
ParKeyVal  *ParKeyVal_add(ParKeyVal *node, const char *key, double *vptr);
double     *ParKeyVal_get(ParKeyVal *node, const char *key);
void        ParKeyVal_print(ParKeyVal *self, FILE *fp);
void        ParKeyVal_sanityCheck(ParKeyVal *self, const char *file, int line);
int         ParKeyVal_equals(ParKeyVal *lhs, ParKeyVal *rhs);
int         legalName(const char *name);

#endif
