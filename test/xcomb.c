/**
   @file xcomb.c
   @brief Test comb.c.

   @copyright Copyright (c) 2019, Alan R. Rogers
   <rogers@anthro.utah.edu>. This file is released under the Internet
   Systems Consortium License, which can be found in file "LICENSE".
*/
#include "comb.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

void usage(void);
int visit(int t, int *c, void *data);
int nvisit(int k, int n[k], int *b[k], void *data);

typedef struct Cdata Cdata;

struct Cdata {
    int verbose;
    int ntot;
    int count;
};

void usage(void) {
    fprintf(stderr,"usage: xcomb [options]\n");
    fprintf(stderr,"where options may include\n");
    fprintf(stderr," -n <a,b,...,z> : number of balls in each box\n");
    fprintf(stderr," -v          : verbose output\n\n");
    fprintf(stderr," Comma-separated arguments to -n must be positive.\n");
    exit(EXIT_FAILURE);
}

int visit(int t, int *c, void *data) {
    Cdata *cdata = (Cdata *) data;
    int n = cdata->ntot, next = 0;
    int i, j;

    cdata->count += 1;

    // print c
    if(cdata->verbose) {
        for(i=0; i<t; ++i) {
            printf("%d", c[i]);
            if(i < t-1)
                putchar(' ');
        }

        if(n > t) {
            // set b equal to complement of combination in c
            j = 0;
            int b[n - t];
            for(i=0; i<t; ++i) {
                while(next < c[i])
                    b[j++] = next++;
                next = c[i] + 1;
            }

            while(j < n - t)
                b[j++] = next++;
            // print b
            fputs(" : ", stdout);
            for(i=0; i < n-t; ++i) {
                printf("%d", b[i]);
                if(i < n-1)
                    putchar(' ');
            }
        }

        putchar('\n');
    }

    return 0;
}

int nvisit(int k, int n[k], int *b[k], void *data) {
    assert(data != NULL);
    Cdata *cdata = (Cdata *) data;
    cdata->count += 1;
    int h, i, j;

    if(cdata->verbose) {
        for(i=0; i<k; ++i) {
            for(j=0; j<n[i]; ++j)
                printf(" %d", b[i][j]);
            fputs(i == k-1 ? "\n" : " :", stdout);
        }
    }

    // total sample size
    int ntot=0;
    for(i=0; i<k; ++i)
        ntot += n[i];

    // Make an 1D array containing all the entries in b.
    int ndx[ntot];
    j=0;
    for(h=0; h<k; ++h)
        for(i=0; i<n[h]; ++i)
            ndx[j++] = b[h][i];

    // sort the array
    qsort(ndx, ntot, sizeof(ndx[0]), compareInts);

    // ndx should now look like [0, 1, 2, ..., ntot-1]
    for(i=0; i<ntot; ++i)
        assert(ndx[i] == i);

    return 0;
}

int main(int argc, char **argv) {
    int i, j, verbose=0, ntot=0;
    int *n = NULL;
    int nboxes = 0;
    char buff[100];
    int  empty_box = 0;

    for(i=1; i<argc; ++i) {
        if(strcmp(argv[i], "-v") == 0)
            verbose=1;
        else if(strcmp(argv[i], "-n") == 0) {
            // parse next arg, which has form "3,1,5"
            i += 1;
            if(i == argc)
                usage();

            // put arg into buff, so it can be modified safely
            strcpy(buff, argv[i]);
            if(*buff == '\0')
                usage();

            // Each comma-separated value is the size of a box.
            // Count the number of boxes.
            nboxes=0;
            for(char *s=buff; *s != '\0'; ++s)
                if(*s == ',')
                    ++nboxes;
            ++nboxes;

            // Allocate an array large enough for the boxes.
            n = malloc(nboxes * sizeof(n[0]));
            CHECKMEM(n);

            // Step through buff, assigning values to n.
            char *token, *next = buff;
            for(j=0; j < nboxes; ++j) {
                token = strsep(&next, ",");
                assert(token != NULL);
                if(*token == '\0') // missing value
                    usage();
                n[j] = strtol(token, NULL, 10);
                if(n[j] < 0)
                    usage();
                if(n[j] == 0)
                    empty_box = 1;
                ntot += n[j];
            }
        }else
            usage();
    }

    // defaults for nboxes, n, and ntot
    if(nboxes == 0) {
        nboxes = 3;
        n = malloc(nboxes * sizeof(n[0]));
        CHECKMEM(n);
        for(i=0; i<nboxes; ++i) {
            n[i] = 2;
            ntot += n[i];
        }
    }

    int x = n[0];

    if(verbose) {
        printf("n:");
        for(i=0; i<nboxes; ++i)
            printf(" %d", n[i]);
        putchar('\n');
        printf("ntot=%d x=%d\n", ntot, n[0]);
    }

    Cdata cdata = {.verbose=verbose, .ntot=ntot, .count=0};

    int status = traverseComb(ntot, x, visit, &cdata);

    int y[2];
    y[0] = n[0];
    y[1] = 0;
    for(i=1; i<nboxes; ++i)
        y[1] += n[i];

    if(verbose) {
        printf("traverseComb returned %d. count=%d\n",
               status, cdata.count);
        printf("Expecting count=%ld\n", binom(ntot, n[0]));
    }

    assert(cdata.count == binom(ntot, n[0]));

    unitTstResult("traverseComb", status==0 ? "OK" : "FAIL");

    cdata.count = 0;

    if(empty_box) {
        printf("Can't test traverseMultiComb because cmd line specified"
               " empty boxes.\n");
        printf("Box sizes:");
        for(i=0; i<nboxes; ++i)
            printf(" %d", n[i]);
        putchar('\n');
        exit(EXIT_FAILURE);
    }

    status = traverseMultiComb(nboxes, n, nvisit, &cdata);

    if(verbose) {
        printf("traverseMultiComb returned %d. count=%d\n",
               status, cdata.count);
        printf("Expecting count=%ld\n", multinom(nboxes, n));
    }

    assert(cdata.count == multinom(nboxes, n));
    unitTstResult("multinom", status==0 ? "OK" : "FAIL");

    assert(binom(0,1) == 0);
    assert(binom(3,4) == 0);
    assert(binom(3,5) == 0);
    assert(binom(5,0) == 1);
    assert(binom(5,1) == 5);
    assert(binom(5,2) == 10);
    assert(binom(5,3) == 10);
    assert(binom(5,4) == 5);
    assert(binom(5,5) == 1);
    assert(binom(0,0) == 1);
    unitTstResult("binom", "OK");

    unitTstResult("traverseMultiComb", status==0 ? "OK" : "FAIL");

    return 0;
}
