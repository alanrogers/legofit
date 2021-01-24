#ifndef ARR_SETPART_H
#define ARR_SETPART_H

#include "typedefs.h"
#include <stdio.h>
#include <stdint.h>

uint64_t      stirling2(uint32_t n, uint32_t k);
void          stirling2_free(void);
long double   lnCoalConst(unsigned n, unsigned k);
double        probPartition(unsigned k, unsigned y[k], long double lnconst);

int traverseSetPartitions(unsigned nelements, unsigned nparts,
                          int (*visit)(unsigned n, unsigned a[n],
                                       void *data), void *data);

#endif
