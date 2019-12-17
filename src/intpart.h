#ifndef ARR_INTPART_H
#define ARR_INTPART_H

int traverseIntPartitions(int n, int m,
                          int (*visit)(int m, int a[m], void *data),
                          void *data);
#endif
