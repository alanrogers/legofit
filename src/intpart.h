#ifndef ARR_INTPART_H
#define ARR_INTPART_H

#include <stdio.h>
#include <stdint.h>

typedef struct NumIntPart NumIntPart;

uint64_t    numIntPart(int32_t n, int32_t k);
void        numIntPart_free(void);
int         traverseIntPartitions(int n, int m,
                                  int (*visit)(int mm, int a[mm], void *data),
                                  void *data);
#endif
