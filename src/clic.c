/**
@file clic.c
@page clic
@author Daniel R. Tabin and Alan R. Rogers
@brief Composite likelihood information criterion.

# `clic`: calculate the composite likelihood information criterion.

This program, like it's sibling @ref bepe "bepe", provides a tool for
selecting among models that differ in complexity. It implements the
"composite likelihood information criterion" (CLIC) of Varin and Vidoni
(2005, Biometrika 92(3):519-528).

CLIC is a generalization of Akaike's information criterion (AIC). It
reduces to -AIC/2 when full likelihood is available. The optimal model
is the one with the lowest value of CLIC.

@copyright Copyright (c) 2018, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "hessian.h"
#include "misc.h"
#include "strdblqueue.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

void usage(void);

const char *usageMsg =
    "usage: clic [options] <f1.pts> <f2.pts> ... -L <f1.legofit> <f2.legofit ...\n"
    " where the .pts files and the .legofit files were all generated by the\n"
    " same set of legofit runs. These runs may be on real data plus bootstrap\n"
    " replicates, or they may be simulation replicates. All files should have\n"
    " the same set of parameters.\n\n"
    "Options:\n"
    "   -h or --help   : print this message\n";

void usage(void) {
    fputs(usageMsg, stderr);
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv){
    int i, j, k;
    int verbose=0, nfiles=0, nLegoFiles=0;
    int gotDashL = 0;
    time_t currtime = time(NULL);

    hdr("clic: composite likelihood information criterion");
#if defined(__DATE__) && defined(__TIME__)
    printf("# Program was compiled: %s %s\n", __DATE__, __TIME__);
#endif
    printf("# Program was run: %s", ctime(&currtime));
    putchar('\n');

    // first pass through arguments counts file names
    for(i=1; i<argc; ++i) {
        if(argv[i][0] == '-') {
            if(strcmp(argv[i], "-L") == 0) {
                gotDashL = 1;
                continue;
            }else {
                fprintf(stderr,"unknown flag argument: %s\n", argv[i]);
                usage();
            }
        }
        if(gotDashL)
            ++nLegoFiles;
        else
            ++nfiles;
    }
    if(nfiles != nLegoFiles) {
        fprintf(stderr, "%s:%d\n"
                " Inconsistent number of files!"
                " %d .pts files and %d legofit files\n",
                __FILE__, __LINE__, nfiles, nLegoFiles);
        usage();
    }
    if(nfiles < 3) {
        fprintf(stderr,"nfiles=%d; need at least 3\n", nfiles);
        usage();
    }

    // 2nd pass builds array of bootstrap filenames
    const char *ptsfname[nfiles], *legofname[nfiles];
    gotDashL=0;
    j=k=0;
    for(i = 1; i < argc; i++) {
        if(argv[i][0] == '-') {
            if(strcmp(argv[i], "-L") == 0) {
                gotDashL = 1;
                assert(k==nfiles);
            }
            continue;
        }
        if(gotDashL)
            legofname[j++] = argv[i];
        else
            ptsfname[k++] = argv[i];
    }
    assert(j==nfiles);
    assert(k==nfiles);

    // Read legofit files into an array of FIFO queues
    StrDblQueue *queue[nfiles];
    for(i=0; i < nfiles; ++i) {
        queue[i] = StrDblQueue_parseLegofit(legofname[i], 1);
        if(queue[i] == NULL) {
            fprintf(stderr,"%s:%d: could not parse legofit file \"%s\"\n",
                    __FILE__,__LINE__, legofname[i]);
            exit(EXIT_FAILURE);
        }
        if(i>0) {
            if(StrDblQueue_compare(queue[0], queue[i])) {
                fprintf(stderr, "%s:%d: inconsistent parameters in"
                        " files %s and %s\n", __FILE__,__LINE__,
                        legofname[0], legofname[i]);
                exit(EXIT_FAILURE);
            }
        }
    }

    // Use the queues to populate an array of parameter names
    // and a matrix of parameter values. Rows are legofit files;
    // columns are parameters.
    const int npar = StrDblQueue_length(queue[0]);
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

    //Check to see if matrix is negative
    if(gsl_matrix_isnonneg (c_matrix))
        fprintf(stderr, "%s: %d: Warning!  Covar matrix isn't negative.\n",
                __FILE__,__LINE__);

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

    // loop over pts files
    printf("#%14s %20s %15s\n", "clic", "PtsFile", "LegofitFile");
    for(k=0; k < nfiles; ++k) {
        // Build hessian matrix from ptsfname
        Hessian hesobj = hessian(ptsfname[k]);
        gsl_matrix *H = hesobj.hessian;
        char **Hparname = hesobj.parname;
        double lnL = hesobj.lnL;

        // Do the two sets of parameters match? (One from legofit files,
        // the other from the pts file.)
        int mismatch = 0;
        if(hesobj.npar != npar) {
            fprintf(stderr,"%s:%d: mismatch in number of parameters\n",
                    __FILE__,__LINE__);
            fprintf(stderr,"    : pts file has %d, legofit files have %d\n",
                    hesobj.npar, npar);
            mismatch=1;
        }
        for(i=0; mismatch==0 && i<npar; ++i) {
            if(0 != strcmp(parname[i], Hparname[i])) {
                fprintf(stderr,"%s:%d: mismatch in parameter names\n",
                        __FILE__,__LINE__);
                mismatch=1;
            }
        }

        if(mismatch) {
            fprintf(stderr,"%s:%d: mismatch between parameters"
                    " in files %s and %s\n", __FILE__,__LINE__,
                    ptsfname[k], legofname[0]);
            exit(EXIT_FAILURE);
        }

        // allocate matrix to hold product H*c_matrix
        gsl_matrix *HC = gsl_matrix_alloc(npar, npar);

        // form matrix product HC = H*c_matrix
        gsl_blas_dgemm(CblasNoTrans, // don't transpose H
                       CblasNoTrans, // don't transpose C
                       1.0,          // scalar
                       H,            // Hermitian matrix
                       c_matrix,     // covariance matrix
                       0.0,          // scalar
                       HC            // product
                       );

        double trace=0.0;
        for(i=0; i<npar; ++i)
            trace += gsl_matrix_get(HC, i, i);

        gsl_matrix_free(H);
        gsl_matrix_free(HC);
        for(i=0; i<npar; ++i)
            free(Hparname[i]);
        free(Hparname);

        /*
          This is the negative of the information criterion of Varin,
          Cristiano and Vidoni, Paolo. 2005. A note on composite
          likelihood inference and model selection. Biometrika
          92(3):519-528. Eqn 5, p 523. We take the negative so that
          the selected model is the one minimizing clic.

          Note that "trace" should be negative at a local maximum. Under
          full likelihood, "trace" is -k, where k is the number of
          parameters, and clic reduces to AIC/2.
        */
        double clic = -(lnL + trace);

        printf("%15.10lg %20s %15s\n", clic,
               mybasename(ptsfname[k]),
               mybasename(legofname[k]));
    }
    
    return 0;
 }
