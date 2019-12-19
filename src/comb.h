#ifndef ARR_COMB_H
#define ARR_COMB_H

int traverseCombinations(int n, int t,
                         int (*visit)(int tt, int a[tt], void *data),
                         void *data);

#endif
