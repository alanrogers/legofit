#include "setpart.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

static int f(int mu, int nu, int sigma, unsigned n, 
              unsigned a[n+1], int (*visit)(unsigned n, unsigned a[n],
                                            void *data), void *data);
static int b(int mu, int nu, int sigma, unsigned n,
              unsigned a[n+1], int (*visit)(unsigned n, unsigned a[n],
                                            void *data), void *data);

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
              unsigned a[n+1], int (*visit)(unsigned n, unsigned a[n],
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
              unsigned a[n+1], int (*visit)(unsigned n, unsigned a[n],
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
            
