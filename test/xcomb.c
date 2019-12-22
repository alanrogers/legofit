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
int nvisit(int k, int n[k], int *b[k], void *data);

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

int nvisit(int k, int n[k], int *b[k], void *data) {
    int i, j;
    for(i=0; i<k; ++i) {
        for(j=0; j<n[i]; ++j)
            printf(" %d", b[i][j]);
        fputs(i == k-1 ? "\n" : " :", stdout);
    }
    return 0;
}

int main(void) {
    int ntot=8, t=4;

    int status = traverseComb(ntot, t, visit, &ntot);

    printf("traverseComb returned %d\n", status);

    int n[] = {1,1,1};
    int k = sizeof(n)/(sizeof(n[0]));
    ntot = 0;
    for(int i=0; i<k; ++i)
        ntot += n[i];

    status = traverseMultiComb(k, n, nvisit, NULL);
    
    printf("traverseMultiComb returned %d\n", status);

    return 0;
}
