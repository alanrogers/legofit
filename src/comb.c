/**
   @file comb.c
   @brief Combinations
   @author Alan R. Rogers
   @copyright Copyright (c) 2019, Alan R. Rogers
   <rogers@anthro.utah.edu>. This file is released under the Internet
   Systems Consortium License, which can be found in file "LICENSE".
*/
#include "comb.h"
#include "misc.h"
#include "u64i64map.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>

typedef struct MCdat MCdat;

MCdat *MCdat_new(int k, int n[k],
                 int (*visit)(int kk, int nn[kk], int *b[kk], void *data),
                 void *data);
void   MCdat_free(MCdat *self);
int    MCvisit(int t, int *a, void *data);

struct MCdat {
    int k;    // number of boxes
    int box;  // current box: 0..(k-2)
    int *n;   // n[i] is number of balls in box i;
    int **b;  // b[i][j] is the index of the j'th ball in box i,
    int cdim; // dimension of array c
    int *c;   // complement: balls not yet assigned to boxes.

    // Pointer to function that will operate on each combination.
    int (*visit)(int kk, int nn[kk], int *b[kk], void *data);

    void *data; // to be manipulated by visit function.
};

MCdat *MCdat_new(int k, int n[k],
                 int (*visit)(int kk, int nn[kk], int *b[kk], void *data),
                 void *data) {
    int i, ntot = 0;
    for(i=0; i<k; ++i)
        ntot += n[i];
    
    MCdat *self = malloc(sizeof(MCdat));
    CHECKMEM(self);
    memset(self, 0, sizeof(MCdat));

    self->k = k;
    self->box = 0;
    self->cdim = ntot;

    self->n = malloc(k * sizeof(self->n[0]));
    CHECKMEM(self->n);

    self->b = malloc(k * sizeof(self->b[0]));
    CHECKMEM(self->b);

    for(i=0; i<k; ++i) {
        self->n[i] = n[i];
        if(n[i] == 0) {
            self->b[i] = NULL;
        }else {
            self->b[i] = malloc(n[i] * sizeof(self->b[i][0]));
            CHECKMEM(self->b[i]);
            memset(self->b[i], 0, n[i] * sizeof(self->b[i][0]));
        }
    }

    // Initially all indices [0,ntot-1) are in c, meaning they are not
    // yet in boxes.
    self->c = malloc(ntot * sizeof(self->c[0]));
    CHECKMEM(self->c);
    for(i=0; i < ntot; ++i)
        self->c[i] = i;

    self->visit = visit;
    self->data = data;

    return self;
}

void MCdat_free(MCdat *self) {
    for(int i=0; i < self->k; ++i) {
        if(self->b[i])
            free(self->b[i]);
    }
    free(self->n);
    free(self->b);
    free(self->c);
    free(self);
}

int MCvisit(int t, int *a, void *data) {
    MCdat *dat = (MCdat *) data;
    CHECKMEM(dat);

    assert(t == dat->n[dat->box]);
    int n = dat->cdim;  // Number of unallocated balls.

    if(n <= t) {
        fprintf(stderr,"%s:%s:%d: Can't choose %d out of %d\n",
                __FILE__,__func__,__LINE__,
                t, n);
        exit(EXIT_FAILURE);
    }
    int cmpl[n - t];    // FAILS IF n==t
    int status, next=0;

    // set cmpl equal to complement of combination in a
    int i, j=0;
    for(i=0; i<t; ++i) {
        while(next < a[i])
            cmpl[j++] = next++;
        next = a[i] + 1;
    }
    while(j < n - t)
        cmpl[j++] = next++;

    // Copy combination into dat->b[dat->box]. Currently, a[i] is the
    // index into dat->c of the i'th ball in the current combination,
    // and dat->c[a[i]] is the index of this ball into the original
    // set of balls.
    for(i=0; i < t; ++i)
        dat->b[dat->box][i] = dat->c[a[i]];

    if(dat->box == dat->k - 2) {
        assert(n-t == dat->n[dat->k - 1]);
        // Copy complement into b[k-1]. Currently, cmpl[i] is the index
        // within dat->c of the i'th ball in the complementary set---
        // i.e. of the balls not in "a". dat->c[cmpl[i]] is the index of
        // this ball in the original set.
        for(i=0; i < n-t; ++i)
            dat->b[dat->k - 1][i] = dat->c[cmpl[i]];

        // Array dat->b now contains a completely specified
        // combination.  Hand it to the visit function whose pointer
        // is in dat.  This will do whatever processing the user wants
        // to do with this combination.
        status = (*dat->visit)(dat->k, dat->n, dat->b, dat->data);
    }else{
        // The combination is not yet completely specified. Proceed to
        // next box. Call to traverseComb will iterate through all
        // ways of allocating things to subsequent boxes, given the
        // contents of previous boxes.

        // Save state
        int box = dat->box;
        int cdim = dat->cdim;
        int c[cdim];
        memcpy(c, dat->c, cdim*sizeof(c[0]));

        // Move to next box
        dat->box += 1;
        assert(dat->box < dat->k);
        dat->cdim = n-t;
        for(i=0; i < dat->cdim; ++i)
            dat->c[i] = c[cmpl[i]];

        // Traverse remaining boxes
        status = traverseComb(dat->cdim, dat->n[dat->box],
                              MCvisit, dat);

        // Restore state
        dat->box = box;
        dat->cdim = cdim;
        memcpy(dat->c, c, cdim*sizeof(c[0]));
    }

    return status;
}

/**
 * Visit each way of allocating N balls among k boxes, such that there
 * are b[0] balls in the first box, b[1] in the second, and so on up to
 * b[k-1], and where N is the sum of the b[i]. For each combination
 * visited, the "visit" function is called.  
 * @param[in] k The number of boxes.
 * @param[in] n[k] The number of balls to be allocated to the k'th box.
 * @param[in] visit A function to be called for each combination
 * visited. The arguments of this function are kk, the number of
 * boxes; nn, an array whose i'th entry is the number of balls in box
 * i; b, an array of kk arrays. The i'th entry of b is an array of
 * nn[i] non-negative integers, whose j'th entry is the index of the
 * j'th ball in box i. These indices range from 0 to N-1 inclusive,
 * and each of these N index values is present exactly once in the
 * two-dimensional array b. The last argument of the visit function is
 * a void pointer, which will hold the address of the data argument
 * passed to the traverseMultiComb function.
 * @param[inout] data If non-NULL, data points to a structure to be
 * manipulated by the visit function.
 */
int traverseMultiComb(int k, int n[k],
                      int (*visit)(int kk, int nn[kk],
                                   int *b[kk], void *data),
                      void *data) {
    if(k < 2) {
        fprintf(stderr,"%s:%s:%d: dimension=%d; must be > 1\n",
                __FILE__,__func__,__LINE__,k);
        exit(EXIT_FAILURE);
    }
    int i, status=0;
    int ntot=0;
    for(i=0; i<k; ++i) {
        if(n[i] <= 0) {
            fprintf(stderr,"%s:%s:%d: number of balls in each box"
                    " must be\n positive, but n[%d]=%d.\n",
                    __FILE__,__func__,__LINE__, i, n[i]);
            exit(EXIT_FAILURE);
        }
        ntot += n[i];
    }
    MCdat *dat = MCdat_new(k, n, visit, data);
    CHECKMEM(dat);

    // This will iterate across all ways of choosing n[0] balls
    // from a total of ntot. In other words, it considers all ways
    // of assigning balls to the 0th box. The MCvisit function
    // will then call traverseMultiComb recursively to deal with
    // the other boxes.
    status = traverseComb(ntot, n[0], MCvisit, dat);

    MCdat_free(dat);
    return status;
}

// Visit all subsets of size t that can be drawn from a larger
// set of size n. Algorithm T, page 359 of Knuth, Donald E. 2011. The
// Art of Computer Programming, Volume 4A.
int traverseComb(int n, int t,
                         int (*visit)(int tt, int *c, void *data),
                         void *data) {

    int status = 0;
    int c[t+3], j, x;

    // The labels in comments below correspond to those in Knuth's
    // pseudocode. 
    for(j=1; j<=t; ++j)  // T1
        c[j] = j-1;
    c[t+1] = n;
    c[t+2] = 0;
    j = t;

    // Knuth's code fails when t == n or t == 0. This corrects the problem.
    if(t == n)
        return (*visit)(t, c+1, data);
    if(t == 0)
        return (*visit)(0, NULL, data);
    
    while(1) {
        if( (status = (*visit)(t, c+1, data)) != 0) // T2
            return status;
        if( j > 0 ) {
            x = j;
        }else{
            if(c[1] + 1 < c[2]) { // T3
                c[1] += 1;
                continue;
            }
            j = 2;
            while(1) { // T4
                assert(j-1 >= 0);
                c[j-1] = j-2;
                assert(j>=0);
                x = c[j] + 1;
                if(x != c[j+1])
                    break;
                j += 1;
            }
            if(j > t)
                break;  // T5: terminate
        }
        c[j] = x; // T6
        j -= 1;
    }

    return status;
}

/**
 * Return N!/(prod x[i]!), the number of ways of allocating N =
 * sum(x[i]) balls among k boxes, with x[i] balls in the ith box.
 *
 * @param[in] k The number of boxes.
 *
 * @param[in] x An array of k positive integers, with x[i]
 * representing the number of balls in box i.
 */
long multinom(int k, int x[k]) {
    int i, n=0;
    for(i=0; i<k; ++i) {
        assert(x[i] > 0);
        n += x[i];
    }
    double ans = lgamma( (double) (n+1));
    for(i=0; i<k; ++i)
        ans -= lgamma( (double) (x[i]+1));

    return (long) floor(exp(ans) + 0.5);
}

// So we don't have to calculate the same value more than once.
static pthread_mutex_t map_lock = PTHREAD_MUTEX_INITIALIZER;
static U64I64Map *map=NULL;

/// Binomial coefficient.
int64_t binom(int32_t n, int32_t x) {
    int status, lockstat;
    uint64_t key;
    int64_t value;

    if(x == 0 || n == x)
        return 1LL;

    if(n == 0) // n==0 && x != 0
        return 0LL;

    if(x < 0)
        return 0LL;

    // Construct a 64-bit key from two 32-bit arguments.
    key = (uint32_t) n;
    key <<= 32;
    key |= (uint32_t) x;

    if(map == NULL) {
        // allocate hash table on first call
        lockstat = pthread_mutex_lock(&map_lock);
        if(lockstat)
            ERR(lockstat, "lock");

        map = U64I64Map_new(128);

        lockstat = pthread_mutex_unlock(&map_lock);
        if(lockstat)
            ERR(lockstat, "unlock");

        CHECKMEM(map);
    }else{
        lockstat = pthread_mutex_lock(&map_lock);
        if(lockstat)
            ERR(lockstat, "lock");

        value = U64I64Map_get(map, key, &status);

        lockstat = pthread_mutex_unlock(&map_lock);
        if(lockstat)
            ERR(lockstat, "unlock");

        if(status == 0)
            return value;
    }

    if(n > 0) {
        value = binom(n-1, x-1) + binom(n-1, x);
    }else{
        long double v = 1.0;
        while(x > 0) {
            v *= n / (long double) x;
            --n;
            --x;
        }
        value = (int64_t) floorl(v + 0.5);
    }

    lockstat = pthread_mutex_lock(&map_lock);
    if(lockstat)
        ERR(lockstat, "lock");

    status = U64I64Map_insert(map, key, value);

    lockstat = pthread_mutex_unlock(&map_lock);
    if(lockstat)
        ERR(lockstat, "unlock");

    if(status) {
        fprintf(stderr,"%s:%d: inserted duplicated value\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    return value;
}

/// Free the hash map used to store binom values.
void binom_free(void) {
    int status = pthread_mutex_lock(&map_lock);
    if(status)
        ERR(status, "lock");

    if(map) {
        U64I64Map_free(map);
        map = NULL;
    }

    status = pthread_mutex_unlock(&map_lock);
    if(status)
        ERR(status, "unlock");
}
