/**
 * @file xstirling2.c
 * @author Alan R. Rogers
 * @brief Test Stirling numbers of the second kind.
 * @copyright Copyright (c) 2019, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "stirling2.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

void usage(void);

void usage(void) {
    fprintf(stderr,"usage: xstirling2 [options]\n");
    fprintf(stderr,"where options may include\n");
    fprintf(stderr," -n <num_samples> : set number of samples\n");
    fprintf(stderr," -v               : verbose output\n");
    fprintf(stderr,"\n   num_samples must be > 1\n");
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {

    int i, verbose = 0;
    long unsigned nmax = 10, n; // overflows at nmax = 27

    for(i=1; i<argc; ++i) {
        if(strcmp(argv[i], "-v") == 0)
            verbose=1;
        else if(strcmp(argv[i], "-n") == 0) {
            i += 1;
            if(i == argc)
                usage();
            nmax = strtol(argv[i], NULL, 10);
        }else
            usage();
    }

    Stirling2 *s = Stirling2_new(nmax);
    CHECKMEM(s);

    if(verbose)
        Stirling2_print(s, stdout);

    assert(Stirling2_val(s, 0, 0) == 1);
    for(n=1; n <= nmax; ++n) {
        assert(Stirling2_val(s,n,0) == 0);
        assert(Stirling2_val(s,n,1) == 1);
        assert(Stirling2_val(s,n,n) == 1);
    }

    if(nmax >= 3)
        assert(Stirling2_val(s, 3, 2) == 3);

    if(nmax >= 4) {
        assert(Stirling2_val(s, 4, 2) == 7);
        assert(Stirling2_val(s, 4, 3) == 6);
    }

    if(nmax >= 5) {
        assert(Stirling2_val(s, 5, 2) == 15);
        assert(Stirling2_val(s, 5, 3) == 25);
        assert(Stirling2_val(s, 5, 4) == 10);
    }

    if(nmax >= 6) {
        assert(Stirling2_val(s, 6, 2) == 31);
        assert(Stirling2_val(s, 6, 3) == 90);
        assert(Stirling2_val(s, 6, 4) == 65);
        assert(Stirling2_val(s, 6, 5) == 15);
    }

    Stirling2_free(s);

    unitTstResult("Stirling2", "OK");
    return 0;
}
