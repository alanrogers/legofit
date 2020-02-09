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
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

typedef struct MCdat MCdat;

MCdat *MCdat_new(int k, int n[k],
                 int (*visit)(int kk, int nn[kk], int *b[kk], void *data),
                 void *data);
void   MCdat_free(MCdat *self);
int    MCvisit(int t, int a[t], void *data);

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
        self->b[i] = malloc(n[i] * sizeof(self->b[i][0]));
        CHECKMEM(self->b[i]);
        memset(self->b[i], 0, n[i] * sizeof(self->b[i][0]));
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
    for(int i=0; i < self->k; ++i)
        free(self->b[i]);
    free(self->n);
    free(self->b);
    free(self->c);
    free(self);
}

int MCvisit(int t, int a[t], void *data) {
    MCdat *dat = (MCdat *) data;
    CHECKMEM(dat);

    assert(t == dat->n[dat->box]);
    int n = dat->cdim;  // Number of unallocated balls.
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
 * a NULL pointer, which will hold the address of the data argument
 * passed to the traverseMultiComb function.
 * @param[inout] data If non-NULL, data points to a structure to be
 * manipulated by the visit function.
 */
int traverseMultiComb(int k, int n[k],
                      int (*visit)(int kk, int nn[kk],
                                   int *b[kk], void *data),
                      void *data) {
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
                         int (*visit)(int tt, int c[tt], void *data),
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

    // Knuth's code fails when t == n. This corrects the problem.
    if(t == n)
        return (*visit)(t, c+1, data);
    
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
    fprintf(stderr,"multinom(%d,x); x=(", k);
    for(i=0; i<k; ++i)
        printf("%d ", x[i]);
    printf(")\n");
    for(i=0; i<k; ++i) {
        assert(x[i] > 0);
        n += x[i];
    }
    double ans = lgamma( (double) (n+1));
    for(i=0; i<k; ++i)
        ans -= lgamma( (double) (x[i]+1));

    return (long) floor(exp(ans) + 0.5);
}

/**
 * Binomial coefficient. Return n!/(x! * (n-x)!), the number of ways
 * of choosing x items out of n.
 *
 * @param[in] n total number of items
 *
 * @param[in] x number to choose out of n.
 */
long binom(long n, long x) {
    double ans = lgamma( (double) (n+1));
    ans -= lgamma( (double) (x+1));
    ans -= lgamma( (double) (n-x+1));

    return (long) floor(exp(ans) + 0.5);
}
