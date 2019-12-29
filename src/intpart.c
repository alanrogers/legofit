/**
   @file intpart.c
   @brief Partitions of an integer into a given number of summands.

   @copyright Copyright (c) 2019, Alan R. Rogers
   <rogers@anthro.utah.edu>. This file is released under the Internet
   Systems Consortium License, which can be found in file "LICENSE".
*/
#include "intpart.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// The number of ways to partition a non-negative integer, n, into k
// positive summands.
struct NumIntPart {
    unsigned nmax;

    // offset is an array of dimension nmax+1, whose i'th entry
    // is offset[i] = (i*(i+1))/2.
    unsigned *offset;

    unsigned nElements; // number of elements in self->s

    // s is a lower triangular matrix (including diagonal), stored as
    // an array. The (n,k)th element of s, where n <= nmax and k <= n,
    // is s[offset[n] + k] and equals P(n,k), the number of ways of
    // partitioning a positive integer, n, into k positive summands.
    unsigned *s;
};

NumIntPart *NumIntPart_new(unsigned nmax) {
    NumIntPart *self = malloc(sizeof(NumIntPart));
    CHECKMEM(self);
    memset(self, 0, sizeof(NumIntPart));

    self->nmax = nmax;

    self->offset = malloc( (nmax+1) * sizeof(self->offset[0]));
    CHECKMEM(self->offset);

    self->nElements = ((nmax+1)*(nmax+2))/2;

    self->s = malloc( (self->nElements) * sizeof(self->s[0]));
    CHECKMEM(self->s);

    unsigned n, k, i, j, ndx;;

    for(n=0; n<=nmax; ++n)
        self->offset[n] = (n*(n+1))/2;

    self->s[0] = 1;
    for(n=1; n<=nmax; ++n)
        self->s[self->offset[n] + 0] = 0;

    for(n=1; n<=nmax; ++n) {
        for(k=1; k<=n; ++k) {
            i = (n-k < k ? 0 : self->s[self->offset[n-k] + k]);
            j = self->s[self->offset[n-1] + k-1];
            ndx = self->offset[n] + k;
            assert(ndx < self->nElements);
            self->s[ndx] = i + j;
        }
    }

    return self;
}

void NumIntPart_free(NumIntPart *self) {
    free(self->offset);
    free(self->s);
    free(self);
}

/// Return the number of ways of partitioning a positive integer n
/// into k summands.
unsigned NumIntPart_val(NumIntPart *self, unsigned n, unsigned k) {
    assert(n <= self->nmax);
    if(!(k <= n)) {
        dostacktrace(__FILE__,__LINE__,stderr);
        fprintf(stderr,"%s:%s:%d: n=%u k=%u\n",
                __FILE__,__func__,__LINE__,n, k);
    }
    assert(k <= n);
    unsigned ndx = self->offset[n] + k;
    assert( ndx < self->nElements );
    return self->s[ndx];
}

void NumIntPart_print(NumIntPart *self, FILE *fp) {
    unsigned i, j;
    fprintf(fp,"%3s:", "n\\k");
    for(j=0; j <= self->nmax; ++j) {
        fprintf(fp, "%5u", j);
        if(j < self->nmax)
            putc(' ', fp);
    }
    putc('\n', fp);
        
    for(i=0; i <= self->nmax; ++i) {
        fprintf(fp, "%3u:", i);
        for(j=0; j <= i; ++j) {
            fprintf(fp, "%5u", NumIntPart_val(self, i, j));
            if(j < self->nmax)
                putc(' ', fp);
        }
        putc('\n', fp);
    }
}

/// Partition a positive integer n into a sum of m positive integers.
/// Algorithm H, p 392 of Knuth, Donald E. 2011. The art of computer
/// programming, volume 4A.
int traverseIntPartitions(int n, int m,
                          int (*visit)(int mm, int a[mm], void *data),
                          void *data) {
    if(m < 2) {
        fprintf(stderr,"%s:%d: m (%d) must be > 2\n",
                __FILE__,__LINE__, m);
        exit(EXIT_FAILURE);
    }
    if(m > n) {
        fprintf(stderr,"%s:%d: m (%d) must be <= n (%d)\n",
                __FILE__,__LINE__, m, n);
        exit(EXIT_FAILURE);
    }
    int j;
    int a[m+1];
    a[0] = n - m + 1;
    for(j=1; j < m; ++j)
        a[j] = 1;
    a[m] = -1;

    while(1) {
        int status = visit(m, a, data);
        if(status)
            return status;
        if(a[1] < a[0] - 1) {
            a[0] -= 1;
            a[1] += 1;
            continue;
        }

        j = 2;
        int s = a[0] + a[1] - 1;
        while(a[j] >= a[0] - 1) {
            s += a[j];
            j += 1;
        }
        if(j+1>m)
            break;
        int x = a[j] + 1;
        a[j] = x;
        j -= 1;
        while(j>0) {
            a[j] = x;
            s -= x;
            j -= 1;
        }
        a[0] = s;
    }
    return 0;
}

