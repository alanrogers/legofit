/**
   @file xcomb.c
   @brief Test comb.c.

   @copyright Copyright (c) 2019, Alan R. Rogers
   <rogers@anthro.utah.edu>. This file is released under the Internet
   Systems Consortium License, which can be found in file "LICENSE".
*/
#include "comb.h"
#include <stdio.h>
#include <assert.h>

int visit(int t, int c[t], void *data);

int visit(int t, int c[t], void *data) {
    int *n = (int *) data;
    int b[*n - t], next=0;

    // set b equal to complement of combination in c
    int i, j=0;
    for(i=0; i<t; ++i) {
        while(next < c[i])
            b[j++] = next++;
        next = c[i] + 1;
    }
    while(j < *n - t)
        b[j++] = next++;

    // print c
    for(i=0; i<t; ++i) {
        printf("%d", c[i]);
        if(i < t-1)
            putchar(' ');
    }

    // print b
    fputs(" : ", stdout);
    for(i=0; i< (*n-t); ++i) {
        printf("%d", b[i]);
        if(i < *n-1)
            putchar(' ');
    }

    putchar('\n');
    return 0;
}

int main(void) {
    int n=8, t=4;

    int status = traverseComb(n, t, visit, &n);

    printf("traverseComb returned %d\n", status);
    return 0;
}
