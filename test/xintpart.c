/**
 * @file xintpart.c
 * @author Alan R. Rogers
 * @brief Generate all partitions of an integer.
 * @copyright Copyright (c) 2019, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "misc.h"
#include "intpart.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

typedef struct VisitDat VisitDat;

// Data manipulated by visit function.
struct VisitDat {
    int verbose;
    int trueSum;
    unsigned long count;
};

void usage(void);
int visit(int m, int a[m], void *data);
int rCmpInts(const void *void_x, const void *void_y);

int visit(int m, int a[m], void *data) {
    VisitDat *vdat = (VisitDat *) data;
    int s=0;
    vdat->count += 1;
    if(vdat->verbose)
        printf("%2lu: ", vdat->count);
    for(int i=0; i < m; ++i) {
        if(vdat->verbose) {
            printf("%d", a[i]);
            if(i < m-1)
                putchar('+');
        }
        s += a[i];
    }
    if(vdat->verbose)
        printf(" = %d", s);

    if(m>0 && vdat->verbose) {
        int d=0, c[m];
        c[0] = 1;
        for(int i=1; i<m; ++i) {
            if(a[i] == a[i-1])
                ++c[d];
            else
                c[++d] = 1;
        }
        ++d;
        qsort(c, d, sizeof(c[0]), rCmpInts);
        int t=0;
        printf(": ");
        for(int i=0; i<d; ++i) {
            printf("%d", c[i]);
            if(i < d-1)
                putchar('+');
            t += c[i];
        }
        printf(" = %d", t);
    }
    putchar('\n');
    return (s == vdat->trueSum ? 0 : 1);
}

void usage(void) {
    fprintf(stderr,"usage: xintpart [options]\n");
    fprintf(stderr,"where options may include\n");
    fprintf(stderr," -n <x>      : integer to be partitioned\n");
    fprintf(stderr," -k <y>      : number of summands per partition\n");
    fprintf(stderr," -v          : verbose output\n\n");
    fprintf(stderr," <x> and <y> must be >0\n");
    fprintf(stderr," <y> must be <= <x>\n");
    exit(EXIT_FAILURE);
}

/**
 * Compare two ints.
 *
 * Function interprets its arguments as pointers to ints.
 *
 * @param void_x,void_y Pointers to the two integers, cast as pointers
 * to voids.
 * @returns -1, 0, or 1 depending on whether the first arg is <,
 * ==, or > the second.
 */
int rCmpInts(const void *void_x, const void *void_y) {
    const int *x = (const int *) void_x;
    const int *y = (const int *) void_y;

    return (*x < *y) ? 1 : (*x > *y) ? -1 : 0;
}

int main(int argc, char **argv) {
    int i, n = 8, k = 3, verbose=0;

    for(i=1; i<argc; ++i) {
        if(strcmp(argv[i], "-v") == 0)
            verbose=1;
        else if(strcmp(argv[i], "-n") == 0) {
            i += 1;
            if(i == argc)
                usage();
            n = strtoul(argv[i], NULL, 10);
        }else if(strcmp(argv[i], "-k") == 0) {
            i += 1;
            if(i == argc)
                usage();
            k = strtoul(argv[i], NULL, 10);
        }else
            usage();
    }

    if(n <= 0 || k <= 0 || k > n)
        usage();

    NumIntPart *nip = NumIntPart_new(n);
    assert(1 == NumIntPart_val(nip, 0, 0));
    for(i=1; i<=n; ++i) {
        assert(0 == NumIntPart_val(nip, i, 0));
        assert(1 == NumIntPart_val(nip, i, i));
    }
    if(n >= 4)
        assert(2 == NumIntPart_val(nip, 4, 2));
    if(verbose)
        NumIntPart_print(nip, stdout);

    VisitDat vdat = {.verbose = verbose,
                     .trueSum = n,
                     .count=0};

    int status = traverseIntPartitions(n, k, visit, &vdat);

    if(status)
        printf("traverseIntPartitions returned %d\n", status);
    if(verbose)
        printf("Found %lu partitions; expected %u\n",
               vdat.count, NumIntPart_val(nip, n, k));

    assert(vdat.count == NumIntPart_val(nip, n, k));
    NumIntPart_free(nip);

    unitTstResult("NumIntPart", "OK");
    unitTstResult("IntPart", status==0 ? "OK" : "FAIL");
    return 0;
}
