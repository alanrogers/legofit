/**
@file clic.c
@page clic
@author Daniel R. Tabin and Alan R. Rogers
@brief Composite likelihood information criterion.

# `clic`: calculate the composite likelihood information criterion.

This program is experimental, and we don't yet recommend using it for
data analysis. The clic statistic involves the expected value of mixed
second partial derivatives of the likelihood function with respect to
pairs of parameter values. Our code calculates this using an
approximation, which we don't yet trust. Our skepticism arises from
the following test. Our composite likelihood function becomes full
likelihood when used with simulated data in which nucleotide sites are
statistically independent. In this case, clic should reduce to
AIC. Yet in our experiments, it does not. We suspect that our
approximation is to blame and are therefore skeptical of our
implementation of clic.

This program, like it's sibling @ref bepe "bepe", provides a tool for
selecting among models that differ in complexity. It implements the
"composite likelihood information criterion" (CLIC) of Varin and Vidoni
(2005, Biometrika 92(3):519-528).

CLIC is a generalization of Akaike's information criterion (AIC). It
reduces to AIC when full likelihood is available and is used in an
analogous fashion.

Our definition of clic is -2 times the criterion defined by Varin and
Vidoni. We multiply by -2 to make our criterion consistent with
AIC. Because of this change, the best model is the one that minimizes
our version of clic. In the version of Varin and Vidoni, the best
model has the maximal value.

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

const char *usageMsg =
     "Usage: clic [options] <file.pts> <boot1.legofit> <boot2.legofit> ...\n"
     "  where file.pts is the .pts file produced by legofit with the real\n"
     "  data, and each \"boot\" file is the legofit output from one bootstrap\n"
     "  replicate. Must include .pts file and at least 2 boostrap files.\n"
     "Options:\n"
     "  -v or --verbose : verbose output\n"
     "  -h or --help    : print this message\n";

void usage(void) {
    fputs(usageMsg, stderr);
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv){
    int i, j;
    int verbose=0;
    int nfiles=0;

    fprintf(stderr,"This program is experimental and should not be trusted.\n");

    // Command line arguments specify file names
    if(argc < 4)
        usage();

    // first pass through arguments counts file names
    for(i=1; i<argc; ++i) {
        if(0==strcmp("-v", argv[i]) || 0==strcmp("--verbose", argv[i]))
            verbose = 1;
        else if(0==strcmp("-h", argv[i]) || 0==strcmp("--help", argv[i]))
            usage();
        else if(argv[i][0] == '-') {
            fprintf(stderr,"Unknown argument: %s\n", argv[i]);
            usage();
        }else
            ++nfiles;
    }
    if(nfiles < 3)
        usage();

    // 2nd pass builds array of bootstrap filenames
    nfiles -= 1;   // counts number of bootstrap files
    const char *ptsfname = NULL;
    const char *bootfname[nfiles];
    int gotPtsFile=0;
    for(i=1, j=0; i<argc; ++i) {
        if(argv[i][0] == '-')
            continue;
        else if(!gotPtsFile) {
            ptsfname = argv[i];
            gotPtsFile=1;
        }else
            bootfname[j++] = argv[i];
    }

    if(ptsfname == NULL)
        usage();

    // Read bootstrap files into an array of FIFO queues
    StrDblQueue *queue[nfiles];
    for(i=0; i < nfiles; ++i) {
        queue[i] = StrDblQueue_parseLegofit(bootfname[i]);
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

    printf("#This program is experimental and should not be trusted.\n");
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
    if(gsl_matrix_isnonneg (c_matrix)) {
        fprintf(stderr, "%s: %d: Warning!  Covar matrix isn't negative."
           "  This is likely due to error!\n", __FILE__,__LINE__);
    }


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

    Hessian hesobj = hessian(ptsfname);
    gsl_matrix *H = hesobj.hessian;
    char **Hparname = hesobj.parname;
    double lnL = hesobj.lnL;

    // Do the two sets of parameters match? (One from bootstrap files,
    // the other from pts file.)
    int mismatch = 0;
    if(hesobj.npar != npar) {
        fprintf(stderr,"%s:%d: mismatch in number of parameters\n",
                __FILE__,__LINE__);
        fprintf(stderr,"    : pts file has %d, boot files have %d\n",
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
                ptsfname, bootfname[0]);
        exit(EXIT_FAILURE);
    }

    // allocate matrix to hold product H*c_matrix
    gsl_matrix *HC = gsl_matrix_alloc(npar,npar);

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

    /*
      This is proportional to the information criterion of Varin,
      Cristiano and Vidoni, Paolo. 2005. A note on composite
      likelihood inference and model selection. Biometrika
      92(3):519-528. Eqn 5, p 523.

      Note that "trace" should be negative at a local maximum. Under
      full likelihood, "trace" is -k, where k is the number of
      parameters, and clic reduces to aic.

      Varin and Vidoni define their information criterion as lnL +
      trace. We multiply by -2 for consistency with AIC. In the
      version of Varin and Vidoni, the best model is the one that
      maximizes their information criterion. In our version, the best
      model minimizes clic.
     */
    double clic = -2.0*(lnL + trace);

    printf("clic: %0.8lg\n %u parameters\n"
           "trace: %lg\n", clic, npar, trace);

    gsl_matrix_free(H);
    gsl_matrix_free(HC);
    for(i=0; i<npar; ++i)
        free(Hparname[i]);
    free(Hparname);
    return 0;
 }
