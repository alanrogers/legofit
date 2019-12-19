/**
   @file comb.c
   @brief Combinations

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
int traverseCombinations(int n, int t,
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
