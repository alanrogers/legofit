#include "hessian.h"
#include "misc.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_errno.h>

Hessian hessian(const char *fname) {

    FILE *fp = fopen(fname, "r");
    if(fp==NULL) {
        fprintf(stderr,"%s:%d: can't read file \"%s\"\n",
                __FILE__,__LINE__, fname);
        exit(EXIT_FAILURE);
    }

    int i, j, k, nfiles, npar, status;
    status = fscanf(fp, "%d %d", &nfiles, &npar);
    if(status != 2)
        DIE("bad return from fscanf");

    // read names of parameters
    char buff[500], parname[npar][50];
    status = fscanf(fp, "%s", buff);
    if(status != 1)
        DIE("bad fscanf");
    
    if(0 != strcmp(buff, "lnL")) {
        fprintf(stderr,"%s:%d: expecting \"lnL\"; got \"%s\"\n",
                __FILE__,__LINE__, stderr);
        exit(EXIT_FAILURE);
    }
    for(i=0; i<npar; ++i) {
        status = fscanf(fp, "%s", parname[i]);
        if(status != 1)
            DIE("bad fscanf");
    }

    // read data
    gsl_vector *lnL = gsl_vector_alloc(nfiles);
    gsl_matrix *data = gsl_matrix_alloc(nfiles, npar);
    double x, y;
    for(i=0; i < nfiles; ++i) {
        status = fscanf(fp, "%lf", &x);
        if(status != 1)
            DIE("bad fscanf");
        gsl_vector_set(lnL, i, x);
        for(j=0; j<npar; ++j) {
            status = fscanf(fp, "%lf", &x);
            if(status != 1)
                DIE("bad fscanf");
            gsl_matrix_set(data, i, j, x);
        }
    }

    // subtract row zero (the maximum-likelihood estimate) to
    // convert to deviations from estimate.
    double estimate[npar];
    for(j=0; j < npar; ++j)
        estimate[j] = gsl_matrix_get(data, 0, j);
    for(i=0; i < nfiles; ++i) {
        for(j=0; j < npar; ++j) {
            x = gsl_matrix_get(data, i, j);
            gsl_matrix_set(data, i, j, x - estimate[j]);
        }
    }

    /*
      Construct design matrix. The first column is filled with 1's and
      represents the Y intercept. The next npar columns have the data
      values from matrix "data". The remaining npar*(npar+1)/2 columns
      have the squares and cross-products of the data values. These
      can be thought of as an unwrapped upper-triangular matrix,
      stored in column-major order. The entry for pair (j,k) is at
      position 1 + npar + offset[k] + j, where offset[k] = k*(k+1)/2.

       0  1  2  3  4  col k
       0  1  3  6 10  offset = k*(k+1)/2
      00 01 02 03 04
      xx 11 12 13 14
      xx xx 22 23 24
      xx xx xx 33 34
      xx xx xx xx 44

      offset[k] = k*(k+1)/2
      x[j][k] is array[offset[k] + j]
    */
    int offset[npar];
    for(j=0; j<npar; ++j)
        offset[j] = (j*(j+1))/2;
    int nterms = 1 + npar + (npar*(npar+1))/2;
    gsl_matrix *design = gsl_matrix_alloc(nfiles, nterms);
    for(i=0; i < nfiles; ++i) {
        gsl_matrix_set(design, i, 0, 1.0); // Y-intercept
        for(j=0; j < npar; ++j) {
            x = gsl_matrix_get(data, i, j);
            gsl_matrix_set(design, i, j+1, x); // linear terms
            for(k=j; k < npar; ++k) {          // quadratic terms
                y = gsl_matrix_get(data, i, k);
                int col = 1 + npar + offset[k] + j;
                gsl_matrix_set(design, i, col, x*y);
            }
        }
    }

    gsl_vector *beta = gsl_vector_alloc(nterms);

    // Fit linear model
    gsl_matrix *cov = gsl_matrix_alloc(nterms, nterms);
    gsl_multifit_linear_workspace *work =
        gsl_multifit_linear_alloc(nfiles, nterms);
    double chisq;
    status = gsl_multifit_linear(design, lnL, beta, cov, &chisq, work);
    if(status) {
        fprintf(stderr,"%s:%d: error %d: %s\n", __FILE__,__LINE__,
                status, gsl_strerror(status));
        exit(EXIT_FAILURE);
    }

    // Construct Hessian matrix
    gsl_matrix *hessian = gsl_matrix_alloc(npar, npar);
    for(j=0; j<npar; ++j) {
        for(k=j; k<npar; ++k) {
            int col = 1 + npar + offset[k] + j;
            double b;
            b = gsl_vector_get(beta, col);
            if(k == j)
                b *= 2;  // 2nd derivative of b*x^2 is 2*b
            gsl_matrix_set(hessian, j, k, b);
            if(j != k) {
                gsl_matrix_set(hessian, k, j, b);
            }
        }
    }

    // Construct return value
    Hessian rval = {
        .npar = npar,
        .lnL = gsl_vector_get(lnL, 0),
        .hessian = hessian
    };
    rval.parname = malloc(npar * sizeof(char *));
    CHECKMEM(rval.parname);
    for(i=0; i<npar; ++i)
        rval.parname[i] = strdup(parname[i]);

    // free memory
    gsl_multifit_linear_free(work);
    gsl_matrix_free(data);
    gsl_matrix_free(design);
    gsl_matrix_free(cov);
    gsl_vector_free(lnL);
    gsl_vector_free(beta);

    return rval;
}
