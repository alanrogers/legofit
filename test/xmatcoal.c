/**
 * @file xmatcoal.c
 * @author Alan R. Rogers
 * @brief Test matcoal.c.
 * @copyright Copyright (c) 2019, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#ifdef NDEBUG
#  error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include "matcoal.h"
#include "misc.h"
#include "mpfrmatcoal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

long double maxAbsErr(int dim, long double x[dim], long double y[dim]);
void usage(void);

long double maxAbsErr(int dim, long double x[dim], long double y[dim]) {
    long double maxerr=0.0, err;

    for(int i=0; i<dim; ++i) {
        err = fabsl(x[i] - y[i]);
        maxerr = fmaxl(maxerr, err);
    }
    return maxerr;
}

void usage(void) {
    fprintf(stderr,"usage: xmatcoal [options]\n");
    fprintf(stderr,"where options may include\n");
    fprintf(stderr," -n <num_samples> : set number of samples\n");
    fprintf(stderr," -v               : verbose output\n");
    fprintf(stderr," -vv              : more verbose\n");
    fprintf(stderr,"\n   num_samples must be > 1\n");
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {

    int nsamples = 15; // overflows at 36, not at 35
    int i, verbose=0;

    for(i=1; i<argc; ++i) {
        if(strcmp(argv[i], "-v") == 0)
            verbose=1;
        else if(strcmp(argv[i], "-vv") == 0)
            verbose=2;
        else if(strcmp(argv[i], "-n") == 0) {
            i += 1;
            if(i == argc)
                usage();
            nsamples = strtol(argv[i], NULL, 10);
        }else
            usage();
    }

    if(nsamples < 2)
        usage();

    assert(0 == MatCoal_nSamples());
    MatCoal_initExterns(nsamples);
    assert(nsamples == MatCoal_nSamples());
    MpfrMatCoal_initExterns(nsamples);
    if(verbose==2) {
        MatCoal_printAll(stdout);
        MpfrMatCoal_printAll(stdout);
    }

    long double ans[nsamples-1], ans2[nsamples-1], err, maxerr=0.0;
    long double v;
    int dim;

    // For each possible dimension (1..(nsamples-1)) and for durations
    // (v) ranging from 0 to infinity, calculate the state vector
    // (using _project) and the expected length of each coalescent
    // interval (using _ciLen). Each calculation is done both in long
    // double precision using MatCoal and in 256-bit precision using
    // MpfrMatCoal. The error is measured as the largest absolute
    // difference between the low-precision and the high-precision
    // calculation. Maxerr is the largest error encountered.
    for(dim=1; dim < nsamples; ++dim) {
        long double eig[dim];
        for(v=0.0; v<10.0; v += 0.5) {
            MatCoal_eigenvals(dim, eig, v);
            MatCoal_project(dim, ans, eig);
            MpfrMatCoal_project(dim, ans2, v);
            err = maxAbsErr(dim, ans, ans2);
            maxerr = fmax(maxerr, err);

            if(dim==2 && fabsl(v - 1.0) < 0.0001) {
                assert( LDbl_near( eig[0], exp(-1.0)));
                assert( LDbl_near( eig[1], exp(-3.0)));
                assert( LDbl_near( ans[0], 1.5*(exp(-1) - exp(-3))));
                assert( LDbl_near( ans[1], exp(-3)));
            }

            if(dim==3 && verbose>1 && fabsl(v-0.5) < 0.0001) {
                printf("project(%d, %Lf):", dim, v);
                for(i=0; i<dim; ++i)
                    printf(" %Lf", ans[i]);
                putchar('\n');
            }

            MatCoal_ciLen(dim, ans, eig);
            MpfrMatCoal_ciLen(dim, ans2, v);
            err = maxAbsErr(dim, ans, ans2);
            maxerr = fmax(maxerr, err);

            if(dim==2 && fabsl(v - 1.0) < 0.0001) {
                assert( LDbl_near( ans[0],
                                  1.0L - 1.5L*expl(-1) + 0.5L*expl(-3.0L)));
                assert( LDbl_near( ans[1],
                                  1.0L/3.0 - expl(-3)/3.0L));
            }

            if(dim==3 && verbose>1 && fabsl(v-0.5L) < 0.0001) {
                printf("ciLen(%d, %Lf):", dim, v);
                for(i=0; i<dim; ++i)
                    printf(" %Lf", ans[i]);
                putchar('\n');
            }

        }

        v = INFINITY;
        MatCoal_eigenvals(dim, eig, v);

        MatCoal_project(dim, ans, eig);
        MpfrMatCoal_project(dim, ans2, v);
        err = maxAbsErr(dim, ans, ans2);
        maxerr = fmaxl(maxerr, err);

        MatCoal_ciLen(dim, ans, eig);
        MpfrMatCoal_ciLen(dim, ans2, v);
        err = maxAbsErr(dim, ans, ans2);
        maxerr = fmaxl(maxerr, err);
    }

    if(verbose)
        fprintf(stderr,"MaxErr = %Le\n", maxerr);

    MatCoal_freeExterns();
    MpfrMatCoal_freeExterns();

    unitTstResult("MatCoal", "OK");
    return 0;
}
