#include <stdio.h>
#include <stdlib.h>

typedef unsigned tipId_t;

static int tipidcmp(const void *vx, const void *vy);

int main(void) {
    const size_t n = 3;
    tipId_t v[n] = {3, 1, 2};

    qsort(v, n, sizeof(v[0]), tipidcmp);

    for(int i=0; i<n; ++i)
        printf(" %u", v[i]);
    putchar('\n');

    return 0;
}

/// This sorts tipId_t values into numerical order, unlike the
/// more complex comparison function, compare_tipId, which is defined
/// in lblndx.c.
static int tipidcmp(const void *vx, const void *vy) {
    tipId_t const * x = vx;
    tipId_t const * y = vy;
    fprintf(stderr,"%s: *x=%u *y=%u\n", __func__, *x, *y);
    if(*x > *y)
        return 1;
    if(*x < *y)
        return -1;
    return 0;
}

