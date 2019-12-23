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

double maxAbsErr(int dim, double x[dim], double y[dim]);
void usage(void);

double maxAbsErr(int dim, double x[dim], double y[dim]) {
    double maxerr=0.0, err;

    for(int i=0; i<dim; ++i) {
        err = fabs(x[i] - y[i]);
        maxerr = fmax(maxerr, err);
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

    MatCoal_initExterns(nsamples);
    MpfrMatCoal_initExterns(nsamples);
    if(verbose==2) {
        MatCoal_printAll(stdout);
        MpfrMatCoal_printAll(stdout);
    }

    double ans[nsamples-1], ans2[nsamples-1], err, maxerr=0.0;
    double v;
    int dim;

    // For each possible dimension (2..(nsamples-1)) and for durations
    // (v) ranging from 0 to infinity, calculate the state vector
    // (using _project) and the expected length of each coalescent
    // interval (using _ciLen). Each calculation is done both in
    // double precision using MatCoal and in 256-bit precision using
    // MpfrMatCoal. The error is measured as the largest absolute
    // difference between the low-precision and the high-precision
    // calculation. Maxerr is the largest error encountered.
    for(dim=2; dim < nsamples; ++dim) {
        for(v=0.0; v<10.0; v += 0.5) {
            MatCoal_project(dim, ans, v);
            MpfrMatCoal_project(dim, ans2, v);
            err = maxAbsErr(dim, ans, ans2);
            maxerr = fmax(maxerr, err);

            MatCoal_ciLen(dim, ans, v);
            MpfrMatCoal_ciLen(dim, ans2, v);
            err = maxAbsErr(dim, ans, ans2);
            maxerr = fmax(maxerr, err);
        }

        v = INFINITY;

        MatCoal_project(dim, ans, v);
        MpfrMatCoal_project(dim, ans2, v);
        err = maxAbsErr(dim, ans, ans2);
        maxerr = fmax(maxerr, err);

        MatCoal_ciLen(dim, ans, v);
        MpfrMatCoal_ciLen(dim, ans2, v);
        err = maxAbsErr(dim, ans, ans2);
        maxerr = fmax(maxerr, err);
    }

    if(verbose)
        fprintf(stderr,"MaxErr = %le\n", maxerr);

    MatCoal_freeExterns();
    MpfrMatCoal_freeExterns();

    unitTstResult("MatCoal", "OK");
    return 0;
}
