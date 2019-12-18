#ifndef ARR_INTPART_H
#define ARR_INTPART_H

#include <stdio.h>

typedef struct NumIntPart NumIntPart;

int         traverseIntPartitions(int n, int m,
                                  int (*visit)(int mm, int a[mm], void *data),
                                  void *data);
NumIntPart *NumIntPart_new(unsigned nmax);
void        NumIntPart_free(NumIntPart *self);
unsigned    NumIntPart_val(NumIntPart *self, unsigned n, unsigned k);
void        NumIntPart_print(NumIntPart *self, FILE *fp);

#endif
