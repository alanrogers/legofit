#ifndef ARR_STIRLING2_H
#define ARR_STIRLING2_H

#include <stdio.h>

typedef struct Stirling2 Stirling2;

Stirling2    *Stirling2_new(long unsigned nmax);
void          Stirling2_free(Stirling2 *self);
long unsigned Stirling2_val(Stirling2 *self, long unsigned n, long unsigned k);
void          Stirling2_print(Stirling2 *self, FILE *fp);

#endif
