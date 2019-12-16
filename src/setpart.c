#include "setpart.h"
#include "stirling2.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

static void generateSetPartitions(SetPart *self);
static void f(int mu, int nu, int sigma, SetPart *self, unsigned n, 
              unsigned a[n+1]);
static void b(int mu, int nu, int sigma, SetPart *self, unsigned n,
              unsigned a[n+1]);
static void visit(SetPart *self, int n, unsigned a[n+1]);

struct SetPart {
    unsigned n; // size of set
    unsigned m; // number of subsets
    long unsigned k; // number of ways to partition n into k subsets

    // A matrix of dimension npart X nelements.  The ij'th entry
    // specifies the subset to which element j belongs in the i'th
    // partition. It is an index into the n-dimensional array
    // representing the original set. For example, suppose that n=4, m
    // = 2, and m[1] = [0,1,0,1]. This means that the 0th and 1th
    // elements belong to subset 0 and the 1th and 3rd to subset 1.
    //
    // The matrix is stored as a 1-dimensional array. The entry for
    // row i and column j is at ndx[i*n + j], where 0 <= i < k, and 0 <=
    // j < n.
    unsigned *ndx;
};

SetPart *SetPart_new(unsigned nelements, unsigned nsubdivisions,
                     Stirling2 *s) {
    SetPart *self = malloc(sizeof(SetPart));
    CHECKMEM(self);
    memset(self, 0, sizeof(SetPart));

    self->n = nelements;
    self->m = nsubdivisions;
    long unsigned npart = Stirling2_val(s, self->n, self->m);
    self->k = 0;
    self->ndx = malloc(npart * self->n * sizeof(self->ndx[0]));
    CHECKMEM(self->ndx);

    generateSetPartitions(self);
    assert(npart == self->k);
    return self;
}

void SetPart_free(SetPart *self) {
    free( self->ndx );
    free(self);
}

unsigned SetPart_sizeOfSet(SetPart *self) {
    return self->n;
}

unsigned SetPart_nSubsets(SetPart *self) {
    return self->m;
}

long unsigned SetPart_nPartitions(SetPart *self) {
    return self->k;
}

/// Fill array p with the indices of the subsets into which
/// the n array elements are assigned by the i'th partition.
void SetPart_getPartition(const SetPart *self, long unsigned i,
                          unsigned n, unsigned p[n]) {
    assert(n == self->n);
    memcpy(p, self->ndx + i*self->n, self->n * sizeof(self->ndx[0]));
}

static void visit(SetPart *self, int n, unsigned a[n+1]) {
    assert(n == self->n);
    // The first entry of a is used by functions f and b but
    // is not part of the answer, so a+1 is the address of
    // the first of the n entries to be copied.
    for(int i=0; i<n; ++i)
        self->ndx[self->k * self->n + i] = a[i+1];
    //    memcpy(self->ndx + self->k*self->n, a+1, n*sizeof(self->ndx[0]));
    self->k += 1;
}

/// Generate all ways of partitioning a set of n objects into m
/// non-empty subsets.  This code implements Knuth's answer to
/// exercise 17 of section 7.2.1.5, on pp 764-765 of The Art of
/// Computer Programming, Volume 4A, part 1. On p 417, Knuth
/// attributes the algorithm to Frank Ruskey [Lecture Notes in
/// Comp. Sci. 762 (1993), 205-206].
static void generateSetPartitions(SetPart *self) {
    unsigned n = self->n;
    unsigned m = self->m;
    unsigned a[n+1];
    memset(a, 0, (n+1)*sizeof(a[0]));
    for(int j=1; j <= m; ++j)
        a[n-m+j] = j-1;
    f(m, n, 0, self, n, a);
}

// Forward recursion
static void f(int mu, int nu, int sigma, SetPart *self, unsigned n, 
              unsigned a[n+1]) {
    if(mu == 2)
        visit(self, n, a);
    else
        f(mu-1, nu-1, (mu+sigma) % 2, self, n, a);
    if( nu == mu + 1) {
        a[mu] = mu - 1;
        visit(self, n, a);
        while(a[nu] > 0) {
            a[nu] -= 1;
            visit(self, n, a);
        }
    }else if(nu > mu+1){
        if( (mu + sigma) % 2 == 1 )
            a[nu-1] = mu-1;
        else
            a[mu] = mu - 1;
        if( (a[nu] + sigma) % 2 == 1)
            b(mu, nu-1, 0, self, n, a);
        else
            f(mu, nu-1, 0, self, n, a);
        while( a[nu] > 0) {
            a[nu] -= 1;
            if( (a[nu] + sigma) % 2 == 1)
                b(mu, nu-1, 0, self, n, a);
            else
                f(mu, nu-1, 0, self, n, a);
        }
    }
}

// Backward recursion
static void b(int mu, int nu, int sigma, SetPart *self, unsigned n,
              unsigned a[n+1]) {
    if( nu == mu+1 ) {
        while( a[nu] < mu-1) {
            visit(self, n, a);
            a[nu] += 1;
        }
        visit(self, n, a);
        a[mu] = 0;
    }else if( nu > mu+1 ) {
        if( (a[nu] + sigma) % 2 == 1)
            f(mu, nu-1, 0, self, n, a);
        else
            b(mu, nu-1, 0, self, n, a);
        while( a[nu] < mu-1 ) {
            a[nu] += 1;
            if( (a[nu] + sigma) % 2 == 1)
                f(mu, nu-1, 0, self, n, a);
            else
                b(mu, nu-1, 0, self, n, a);
        }
        if( (mu + sigma) % 2 == 1)
            a[nu - 1] = 0;
        else
            a[mu] = 0;
    }
    if( mu == 2)
        visit(self, n, a);
    else
        b(mu-1, nu-1, (mu+sigma) % 2, self, n, a);
}
            
