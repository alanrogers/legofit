#ifndef ARR_SETPART_H
#define ARR_SETPART_H

#include <stdio.h>
#include <stdint.h>

typedef struct Stirling2 Stirling2;

uint64_t      stirling2(uint32_t n, uint32_t k);
long double   lnCoalConst(unsigned n, unsigned k);
double        probPartition(unsigned k, unsigned y[k], long double lnconst);

int traverseSetPartitions(unsigned nelements, unsigned nparts,
                          int (*visit)(unsigned n, unsigned a[n],
                                       void *data), void *data);

#endif
