#ifndef ARR_COMB_H
#define ARR_COMB_H

#include <stdint.h>

int traverseComb(int n, int t,
                 int (*visit)(int tt, int a[tt], void *data),
                 void *data);

int traverseMultiComb(int k, int n[k],
                      int (*visit)(int kk, int nn[kk],
                                   int *b[kk], void *data),
                      void *data);
long multinom(int k, int x[k]);
int64_t binom(int32_t n, int32_t x);
void binom_free(void);

#endif
