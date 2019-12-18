/**
 * @file xsetpart.c
 * @author Alan R. Rogers
 * @brief Generate all partitions of a set of n elements into m subsets.
 * @copyright Copyright (c) 2019, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "setpart.h"
#include "stirling2.h"
#include "misc.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

void usage(void);
int visit(unsigned n, unsigned a[n], void *data);

typedef struct VisitDat VisitDat;

// Data manipulated by visit function.
struct VisitDat {
    int verbose;
    unsigned nsubdivisions;
    unsigned long count;
};

/// Function to process partition. This function prints
/// the partition if verbose!=0, checks that it is valid,
/// and counts the partitions. Returns non-zero if an error is found.
int visit(unsigned n, unsigned a[n], void *data) {
    VisitDat *vdat = (VisitDat *) data;
    if(vdat->verbose) {
        for(int i=0; i<n; ++i) {
            printf("%d", a[i]);
            if(i < n-1)
                putchar(' ');
        }
        putchar('\n');
    }
    unsigned max = 0;
    int status=0;
    for(int i=0; i<n; ++i) {
        if(a[i] > max) {
            if(a[i] != max+1)
                ++status;
            max = a[i];
        }
    }
    if(max != 1 + vdat->nsubdivisions)
        ++status;
    vdat->count += 1;
    return status;
}

void usage(void) {
    fprintf(stderr,"usage: xsetpart [options]\n");
    fprintf(stderr,"where options may include\n");
    fprintf(stderr," -n <n_elements> : size of set\n");
    fprintf(stderr," -m <n_subdiv>   : number of subdivisions\n");
    fprintf(stderr," -v              : verbose output\n\n");
    fprintf(stderr,"   n_elements must be > 0\n");
    fprintf(stderr,"   n_subdiv must be > 0, <= n_elements\n");
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
    int i, verbose = 0;
    unsigned nelem = 6, nsub = 3;

    for(i=1; i<argc; ++i) {
        if(strcmp(argv[i], "-v") == 0)
            verbose=1;
        else if(strcmp(argv[i], "-n") == 0) {
            i += 1;
            if(i == argc)
                usage();
            nelem = strtoul(argv[i], NULL, 10);
        }else if(strcmp(argv[i], "-m") == 0) {
            i += 1;
            if(i == argc)
                usage();
            nsub = strtoul(argv[i], NULL, 10);
        }else
            usage();
    }

    VisitDat vdat = {.verbose = verbose, 
                     .nsubdivisions = nsub,
                     .count=0};

    int status = traverseSetPartitions(nelem, nsub, visit, &vdat);

    Stirling2 *s = Stirling2_new(nelem);
    CHECKMEM(s);
    assert(vdat.count == Stirling2_val(s, nelem, nsub));

    unitTstResult("SetPart", status==0 ? "OK" : "FAIL");
    return 0;
}
