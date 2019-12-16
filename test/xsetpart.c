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

void usage(void);

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
    unsigned nelem = 6, nsub = 3; // overflows at nmax = 27

    for(i=1; i<argc; ++i) {
        if(strcmp(argv[i], "-v") == 0)
            verbose=1;
        else if(strcmp(argv[i], "-n") == 0) {
            i += 1;
            if(i == argc)
                usage();
            nelem = strtol(argv[i], NULL, 10);
        }else if(strcmp(argv[i], "-m") == 0) {
            i += 1;
            if(i == argc)
                usage();
            nsub = strtol(argv[i], NULL, 10);
        }else
            usage();
    }

    Stirling2 *stirling2 = Stirling2_new(nelem);
    CHECKMEM(stirling2);

    SetPart *sp = SetPart_new(nelem, nsub, stirling2);
    CHECKMEM(sp);

    long unsigned k, npart = SetPart_nPartitions(sp);

    assert(nelem == SetPart_sizeOfSet(sp));
    assert(nsub == SetPart_nSubsets(sp));
    assert(npart == Stirling2_val(stirling2, nelem, nsub));

    unsigned p[nelem];

    for(k=0; k < npart; ++k) {
        SetPart_getPartition(sp, k, nelem, p);
        if(verbose) {
            for(i=0; i < nelem; ++i) {
                printf("%d", p[i]);
                if(i < nelem-1)
                    putchar(' ');
            }
            putchar('\n');
        }
        assert(p[0] == 0);
        unsigned max = 0;
        for(i=1; i<nelem; ++i) {
            if(p[i] > max) {
                assert(p[i] = max+1);
                max = p[i];
            }
        }
        assert(max == nsub - 1);
    }
    SetPart_free(sp);
    unitTstResult("SetPart", "OK");
    return 0;
}
