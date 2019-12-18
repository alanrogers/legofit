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
        printf(" = %d\n", s);
    return (s == vdat->trueSum ? 0 : 1);
}

void usage(void) {
    fprintf(stderr,"usage: xintpart [options]\n");
    fprintf(stderr,"where options may include\n");
    fprintf(stderr," -n <x>      : integer to be partitioned\n");
    fprintf(stderr," -m <y>      : number of summands per partition\n");
    fprintf(stderr," -v          : verbose output\n\n");
    fprintf(stderr," <x> and <y> must be >0\n");
    fprintf(stderr," <y> must be <= <x>\n");
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
    int i, n = 8, m = 3, verbose=0;

    for(i=1; i<argc; ++i) {
        if(strcmp(argv[i], "-v") == 0)
            verbose=1;
        else if(strcmp(argv[i], "-n") == 0) {
            i += 1;
            if(i == argc)
                usage();
            n = strtoul(argv[i], NULL, 10);
        }else if(strcmp(argv[i], "-m") == 0) {
            i += 1;
            if(i == argc)
                usage();
            m = strtoul(argv[i], NULL, 10);
        }else
            usage();
    }

    if(n <= 0 || m <= 0 || m > n)
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

    int status = traverseIntPartitions(n, m, visit, &vdat);

    if(status)
        printf("traverseIntPartitions returned %d\n", status);
    if(verbose)
        printf("Found %lu partitions.\n", vdat.count);

    assert(vdat.count == NumIntPart_val(nip, n, m));
    NumIntPart_free(nip);

    unitTstResult("NumIntPart", "OK");
    unitTstResult("IntPart", status==0 ? "OK" : "FAIL");
    return 0;
}
