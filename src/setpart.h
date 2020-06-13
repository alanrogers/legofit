#ifndef ARR_SETPART_H
#define ARR_SETPART_H

#include <stdio.h>

typedef struct Stirling2 Stirling2;

Stirling2    *Stirling2_new(long unsigned nmax);
void          Stirling2_free(Stirling2 *self);
long unsigned Stirling2_val(Stirling2 *self, long unsigned n, long unsigned k);
void          Stirling2_print(Stirling2 *self, FILE *fp);
long double   lnCoalConst(unsigned n, unsigned k);
double        probPartition(unsigned k, unsigned y[k], long double lnconst);

int traverseSetPartitions(unsigned nelements, unsigned nparts,
                          int (*visit)(unsigned n, unsigned a[n],
                                       void *data), void *data);

#endif
