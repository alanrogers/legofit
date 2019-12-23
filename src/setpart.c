#include "setpart.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

static int f(int mu, int nu, int sigma, unsigned n, 
              unsigned a[n+1], int (*visit)(unsigned nn, unsigned a[nn],
                                            void *data), void *data);
static int b(int mu, int nu, int sigma, unsigned n,
              unsigned a[n+1], int (*visit)(unsigned nn, unsigned a[nn],
                                            void *data), void *data);

// Stirling numbers of the second kind.
// S(n,k) is the number of ways of partitioning a set of n objects
// into k subsets.
struct Stirling2 {
    long unsigned nmax;  // overflows at 27

    // offset is an array of dimension nmax+1, whose i'th entry
    // is offset[i] = (i*(i+1))/2.
    long unsigned *offset;

    long unsigned nElements; // number of elements in self->s

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

    self->nElements = ((nmax+1)*(nmax+2))/2;

    self->s = malloc( (self->nElements) * sizeof(self->s[0]));
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
            long unsigned ndx = self->offset[n] + k;
            assert(ndx < self->nElements);
            self->s[ndx] = k * i + j;
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
    if(!(k <= n)) {
        dostacktrace(__FILE__,__LINE__,stderr);
        fprintf(stderr,"%s:%s:%d: n=%lu k=%lu\n",
                __FILE__,__func__,__LINE__,n, k);
    }
    assert(k <= n);
    long unsigned ndx = self->offset[n] + k;
    assert( ndx < self->nElements );
    return self->s[ndx];
}

void Stirling2_print(Stirling2 *self, FILE *fp) {
    long unsigned i, j;
    fprintf(fp,"%3s:", "n\\k");
    for(j=0; j <= self->nmax; ++j) {
        fprintf(fp, "%5lu", j);
        if(j < self->nmax)
            putc(' ', fp);
    }
    putc('\n', fp);
        
    for(i=0; i <= self->nmax; ++i) {
        fprintf(fp, "%3lu:", i);
        for(j=0; j <= i; ++j) {
            fprintf(fp, "%5lu", Stirling2_val(self, i, j));
            if(j < self->nmax)
                putc(' ', fp);
        }
        putc('\n', fp);
    }
}

/**
   Generate all ways of partitioning a set of n objects into m
   non-empty subsets.  This code implements Knuth's answer to
   exercise 17 of section 7.2.1.5, on pp 764-765 of The Art of
   Computer Programming, Volume 4A, part 1. On p 417, Knuth
   attributes the algorithm to Frank Ruskey [Lecture Notes in
   Comp. Sci. 762 (1993), 205-206].

   The "visit" function is called for each partition. Its arguments
   include the current partition and a pointer to a data structure.
   The visit function and the data structure should be defined by the
   user.
*/
int traverseSetPartitions(unsigned nelements, unsigned nsubdivisions,
                           int (*visit)(unsigned n, unsigned a[n],
                                        void *data), void *data) {
    unsigned a[nelements+1];
    memset(a, 0, (nelements+1)*sizeof(a[0]));
    for(int j=1; j <= nsubdivisions; ++j)
        a[nelements-nsubdivisions+j] = j-1;
    if(nsubdivisions==1) {
        (*visit)(nelements, a+1, data);
        return 0;
    }
    return f(nsubdivisions, nelements, 0, nelements, a, visit, data);
}

// Forward recursion
static int f(int mu, int nu, int sigma, unsigned n, 
              unsigned a[n+1], int (*visit)(unsigned nn, unsigned a[nn],
                                            void *data), void *data) {
    int status;
    if(mu == 2)
        (*visit)(n, a+1, data);
    else {
        status = f(mu-1, nu-1, (mu+sigma) % 2, n, a, visit, data);
        if(status)
            return status;
    }
    if( nu == mu + 1) {
        a[mu] = mu - 1;
        (*visit)(n, a+1, data);
        while(a[nu] > 0) {
            a[nu] -= 1;
            (*visit)(n, a+1, data);
        }
    }else if(nu > mu+1){
        if( (mu + sigma) % 2 == 1 )
            a[nu-1] = mu-1;
        else
            a[mu] = mu - 1;
        if( (a[nu] + sigma) % 2 == 1) {
            status = b(mu, nu-1, 0, n, a, visit, data);
            if(status)
                return status;
        }else{
            status = f(mu, nu-1, 0, n, a, visit, data);
            if(status)
                return status;
        }
        while( a[nu] > 0) {
            a[nu] -= 1;
            if( (a[nu] + sigma) % 2 == 1) {
                status = b(mu, nu-1, 0, n, a, visit, data);
                if(status)
                    return status;
            }else{
                status = f(mu, nu-1, 0, n, a, visit, data);
                if(status)
                    return status;
            }
        }
    }
    return 0;
}

// Backward recursion
static int b(int mu, int nu, int sigma, unsigned n,
              unsigned a[n+1], int (*visit)(unsigned nn, unsigned a[nn],
                                            void *data), void *data) {
    int status;
    if( nu == mu+1 ) {
        while( a[nu] < mu-1) {
            (*visit)(n, a+1, data);
            a[nu] += 1;
        }
        (*visit)(n, a+1, data);
        a[mu] = 0;
    }else if( nu > mu+1 ) {
        if( (a[nu] + sigma) % 2 == 1) {
            status = f(mu, nu-1, 0, n, a, visit, data);
            if(status)
                return status;
        }else{
            status = b(mu, nu-1, 0, n, a, visit, data);
            if(status)
                return status;
        }
        while( a[nu] < mu-1 ) {
            a[nu] += 1;
            if( (a[nu] + sigma) % 2 == 1) {
                status = f(mu, nu-1, 0, n, a, visit, data);
                if(status)
                    return status;
            }else{
                status = b(mu, nu-1, 0, n, a, visit, data);
                if(status)
                    return status;
            }
        }
        if( (mu + sigma) % 2 == 1)
            a[nu - 1] = 0;
        else
            a[mu] = 0;
    }
    if( mu == 2)
        (*visit)(n, a+1, data);
    else {
        status = b(mu-1, nu-1, (mu+sigma) % 2, n, a, visit, data);
        if(status)
            return status;
    }
    return 0;
}
            
