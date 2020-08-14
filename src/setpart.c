#include "setpart.h"
#include "misc.h"
#include "u64u64map.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

static int f(int mu, int nu, int sigma, unsigned n, 
              unsigned a[n+1], int (*visit)(unsigned nn, unsigned a[nn],
                                            void *data), void *data);
static int b(int mu, int nu, int sigma, unsigned n,
              unsigned a[n+1], int (*visit)(unsigned nn, unsigned a[nn],
                                            void *data), void *data);

// So we don't have to calculate the same value more than once.
static U64U64Map *map=NULL;

// Stirling numbers of the second kind.  stirling2(n,k) is the number
// of ways of partitioning a set of n objects into k subsets.
uint64_t stirling2(uint32_t n, uint32_t k) {
    uint64_t key, value;
    int status;

    // Initial conditions
    if(n==0 && k==0)
        return 1ULL;
    
    if(n==0 || k==0)
        return 0ULL;

    // Construct a 64-bit key from two 32-bit arguments.
    key = n;
    key <<= 32;
    key |= k;

    if(map == NULL) {
        // allocate hash table on first call
        map = U64U64Map_new();
        CHECKMEM(map);
    }else{
        value = U64U64Map_get(map, key, &status);
        if(status == 0)
            return value;
    }

    // Fundamental recurrence relation
    value = k*stirling2(n-1, k) + stirling2(n-1, k-1);
    
    status = U64U64Map_insert(map, key, value);
    if(status) {
        fprintf(stderr,"%s:%d: inserted duplicated value\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    return value;
}

/// Free the hash map used to store stirling2 values.
void stirling2_free(void) {
    if(map) {
        U64U64Map_free(map);
        map = NULL;
    }
}

/** 
 * Log of constant in coalescent probability from theorem 1.5, p. 11,
 * Durrett, Richard. 2008. Probability Models for DNA Sequence
 * Evolution.
 */
long double lnCoalConst(unsigned n, unsigned k) {
    assert(n >= k);
    assert(k > 0);
    return lgammal(k+1)
        - lgammal(n+1)
        + lgammal(n-k+1)
        + lgammal(k)
        - lgammal(n);
}

/**
 * Probability of a set partition under the coalescent process. There
 * are n gene copies (descendants) at the recent end of a segment of
 * population history. At some point, earlier in the segment, there
 * are k ancestral gene copies. These ancestors define a partition of
 * the set of descendants into k subsets, each representing the
 * descendants of a single ancestor. This function calculates the
 * probability of a given set partition.
 *
 * See theorem 1.5, p. 11, Durrett, Richard. 2008. Probability Models
 * for DNA Sequence Evolution.
 *
 * @param[in] k number of ancestral gene copies
 * @param[in] y[k] array whose i'th entry is the number of descendants
 * of the i'th ancestor. The sum of y[i] is n, the number of
 * descendants. 
 * @param[in] lnconst portion of the probability that does not depend
 * on y and can be calculated in advance, using function lnCoalConst.
 */
double probPartition(unsigned k, unsigned y[k], long double lnconst) {
    assert(k > 0);
    long double x = 0.0L;
    for(unsigned i=0; i<k; ++i)
        x += lgammal(y[i] + 1);
    return (double) expl(lnconst + x);
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
int traverseSetPartitions(unsigned nelements, unsigned nparts,
                           int (*visit)(unsigned n, unsigned a[n],
                                        void *data), void *data) {
    unsigned a[nelements+1];
    memset(a, 0, (nelements+1)*sizeof(a[0]));
    for(int j=1; j <= nparts; ++j)
        a[nelements-nparts+j] = j-1;
    if(nparts==1) {
        (*visit)(nelements, a+1, data);
        return 0;
    }
    return f(nparts, nelements, 0, nelements, a, visit, data);
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

