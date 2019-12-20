/**
   @file comb.c
   @brief Combinations
   @author Alan R. Rogers
   @copyright Copyright (c) 2019, Alan R. Rogers
   <rogers@anthro.utah.edu>. This file is released under the Internet
   Systems Consortium License, which can be found in file "LICENSE".
*/
#include "comb.h"
#include <stdio.h>
#include <assert.h>

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
 * Visit each way of allocating N balls among k boxes, such that there
 * are b[0] balls in the first box, b[1] in the second, and so on up to
 * b[k-1], and where N is the sum of the b[i]. For each combination
 * visited, the "visit" function is called.  
 * @param[in] k The number of boxes.
 * @param[in] b[k] The number of balls to be allocated to the k'th box.
 * @param[in] visit A function to be called for each combination
 * visited. The arguments of this function are kk, the number of
 * boxes; bb, an array whose i'th entry is the number of balls in box
 * i; c, an array of kk arrays. The i'th entry of c is an array of
 * b[i] non-negative integers, whose j'th entry is the index of the
 * j'th ball in box i. These indices range from 0 to N-1 inclusive,
 * and each of these N index values is present exactly once in the
 * two-dimensional array c. The last argument of the visit function is
 * a NULL pointer, which will hold the address of the data argument
 * passed to the traverseMultiComb function.
 * @param[inout] data If non-NULL, data points to a structure to be
 * manipulated by the visit function.
 */
int traverseMultiComb(int k, int b[k],
                      int (*visit)(int kk, int bb[kk],
                                   int *c[kk], void *data),
                      void *data) {
    int i, status=0;
    int n=0;
    for(i=0; i<k; ++i)
        n += b[i];
    return status;
}
