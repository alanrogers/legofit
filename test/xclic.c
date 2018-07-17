/**
@file clic.c
@page clic
@author Daniel R. Tabin
@brief Composite likelihood information criterion.

# `clic`: calculate the composite likelihood information criterion.

This program program will test clic.  More information can be found
in ../src/clic.c

@copyright Copyright (c) 2018, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "hessian.h"
#include "strdblqueue.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

void usage(void);

const char *usageMsg = "usage: xstrdblstck [-v]\n";

void usage(void) {
    fputs(usageMsg, stderr);
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv){

    //setup

    const char *fakeBoot =
    "Fitted parameter values\n"
    "    9 free:\n"
    "        Tnd = 1\n"
    "        Txy = 2\n"
    "         Td = 3\n"
    "         Ta = 4\n"
    "       2Nnd = 5\n"
    "        2Nn = 6\n"
    "       2Nxy = 7\n"
    "     2Nxynd = 8\n"
    "         mN = 9\n"
    "    0 constrained:\n"
    "#       SitePat  BranchLen\n"
    "              x 38177.7670041\n"
    "              y 36976.1129329\n"
    "              n 39773.2585902\n"
    "              d 41503.9037230\n"
    "            x:y 22741.5617898\n"
    "            x:n 2973.2394928\n"
    "            x:d 3177.9284049\n"
    "            y:n 3699.1331968\n"
    "            y:d 3002.9907387\n"
    "            n:d 19125.0587195\n"
    "          x:y:n 7062.0871479\n"
    "          x:y:d 6892.8717880\n"
    "          x:n:d 5889.0132020\n"
    "          y:n:d 6539.7112355\n";

    FILE       *fakeFile = fopen("s1boot0.legofit", "w");
    assert(fakeFile);
    fputs(fakeBoot, fakeFile);
    fclose(fakeFile);

    bool verbose = 0;

      switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0) {
            usage();
        }
        verbose = 1;
        break;
    default:
        usage();
    }

    int files = 2;

    //test matrix

    StrDblQueue *queue[nfiles];
    for(i=0; i < nfiles; ++i) {
        queue[i] = parseLegofit_CLIC("s1boot0.legofit");
        if(i>0) {
            if(StrDblQueue_compare(queue[0], queue[i])) {
                fprintf(stderr, "%s:%d: inconsistent parameters in"
                        " files %s and %s\n", __FILE__,__LINE__,
                        bootfname[0], bootfname[i]);
                exit(EXIT_FAILURE);
            }
        }
    }

    // Use the queues to populate an array of parameter names
    // and a matrix of parameter values. Rows are bootstrap
    // replicates; columns are parameters.
    int npar = StrDblQueue_length(queue[0]);
    char *parname[npar];
    double datmat[nfiles][npar];
    for(i=0; i < nfiles; ++i) {
        for(j=0; j < npar; ++j) {
            StrDbl strdbl;
            queue[i] = StrDblQueue_pop(queue[i], &strdbl);
            datmat[i][j] = strdbl.val;
            if(i==0)
                parname[j] = strdup(strdbl.str);
        }
        assert(queue[i] == NULL); // check that queues are freed
    }

    if(verbose) {
        // Print data matrix with column header
        for(j=0; j < npar; ++j)
            printf(" %s", parname[j]);
        putchar('\n');
        for(i=0; i<nfiles; ++i) {
            for(j=0; j < npar; ++j)
                printf(" %lg", datmat[i][j]);
            putchar('\n');
        }
    }

    // Make covariance matrix
    gsl_matrix *c_matrix = gsl_matrix_alloc(npar,npar);
    make_covar_matrix(nfiles, npar, datmat, c_matrix);

    if(verbose) {
        // Print it
        for (j = 0; j < npar; j++)
            printf(" %8s", parname[j]);
        putchar('\n');
        for (i = 0; i < npar; i++){
            for (j = 0; j < npar; j++){
                printf(" %8.2lg", gsl_matrix_get(c_matrix, i, j));
            }
            printf("\n");
        }
    }

    //test points

    //finish
    remove("s1boot0.legofit");

    printf("All tests completed\n");

 }
