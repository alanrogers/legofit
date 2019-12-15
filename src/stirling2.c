#include "misc.h"
#include <stdio.h>
#include <string.h>

typedef struct Stirling2 Stirling2;

void Stirling2_free(Stirling2 *self);
Stirling2 *Stirling2_new(long unsigned nmax);
long unsigned Stirling2_val(Stirling2 *self, long unsigned n, long unsigned k);

// Stirling numbers of the second kind.
// S(n,k) is the number of ways of partitioning a set of n objects
// into k subsets.
struct Stirling2 {
    long unsigned nmax;  // overflows at 27

    // offset is an array of dimension nmax+1, whose i'th entry
    // is offset[i] = (i*(i+1))/2.
    long unsigned *offset; 

    // s is a lower triangular matrix (including diagonal), stored as
    // an array. The (n,k)th element of s, where n <= nmax and k <= n,
    // is s[offset[n] + k] and equals S(n,k), the number of ways of
    // partitioning a set of n objects into k subsets.
    long unsigned *s;
};

Stirling2 *Stirling2_new(long unsigned nmax) {
    Stirling2 *self = malloc(sizeof(Stirling2));
    CHECKMEM(self);
    memset(self, 0, sizeof(Stirling2));

    self->nmax = nmax;

    self->offset = malloc( (nmax+1) * sizeof(self->offset[0]));
    CHECKMEM(self->offset);

    self->s = malloc( (((nmax+1)*(nmax+2))/2) * sizeof(self->s[0]));
    CHECKMEM(self->s);

    long unsigned n, k;

    for(n=0; n<=nmax; ++n)
        self->offset[n] = (n*(n+1))/2;

    self->s[0] = 1;
    for(n=1; n<=nmax; ++n)
        self->s[self->offset[n] + 0] = 0;

    for(n=1; n<=nmax; ++n) {
        for(k=1; k<=n; ++k) {
            long unsigned i = (k > n-1 ? 0 : self->s[self->offset[n-1] + k]);
            long unsigned j = (k > n ? 0 : self->s[self->offset[n-1] + k-1]);
            self->s[self->offset[n] + k] = k * i + j;
        }
    }

    return self;
}

void Stirling2_free(Stirling2 *self) {
    free(self->offset);
    free(self->s);
    free(self);
}

/// Return S(n,k), the number of ways of partitioning a set
/// of n elements into k subsets.
long unsigned Stirling2_val(Stirling2 *self, long unsigned n, long unsigned k) {
    assert(n <= self->nmax);
    assert(k <= n);
    return self->s[self->offset[n] + k];
}

int main(void) {
    // overflows at nmax = 27
    long unsigned nmax = 10, n, k;

    Stirling2 *s = Stirling2_new(nmax);
    CHECKMEM(s);

    for(n=0; n<=nmax; ++n) {
        for(k=0; k<=n; ++k)
            printf(" %5lu", Stirling2_val(s, n, k));
        putchar('\n');
    }
    return 0;
}
