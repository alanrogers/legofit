#include "rational.h"
#include "matcoal.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>

static void init_dim(int nLin);

// Number of haploid samples in the data. 0 until initialized
static int nsamples = 0;

/*
  Layout of upper triangular matrices is row major:

  00 01 02 03 04  row 0: offset = 0
  xx 11 12 13 14  row 1: offset = dim-1
  xx xx 22 23 24  row 2: offset = 2*dim - 3
  xx xx xx 33 34  row 3: offset = 3*dim - 6
  xx xx xx xx 44  row 4: offset = 4*dim - 10

  offset[i] = i*dim - i*(i+1)/2 = ((2*dim - 1 - i)*i)/2
  x[i][j] is array[offset[i] + j]

  Number of stored elements is dim*(dim+1)/2.
 */

// Array of pointers to matrices of scaled column eigenvectors
// Each matrix is upper triangular and is stored as an array.
// cmat[k] is a matrix of dimension (k+2)X(k+2) and is used
// in epochs with k+3 lineages.
static double **cmat = NULL;

// offset[k] is a pointer to an array of k offsets, which are
// used to address elements of the triangular matrix cmat[k].
static int **offset;

// Array of pointers to matrices for calculating expected branch
// lengths. The kth matrix is rectangular, with dimension (k+2) X
// (k+3), but allocated as a linear array. (i,j)th element of k'th
// matrix is bmat[k][i*(k+3) + j]
static double **bmat = NULL;

// The i'th entry is (i+2) choose 2
static double *beta = NULL;

// Initialize arrays for epochs in which there are nLin Lineages.
// For an epoch with n lineages, dim=n-1, and k=n-2.
static void init_dim(int dim) {
    int i, j, ii, jj;
    int k = dim-1; // index into external arrays
    int nLin = dim+1; // number of lineages in epoch.
    size_t size = sizeof(cmat[k][0]) * (dim*(dim+1))/2;
    cmat[k] = malloc(size);
    CHECKMEM(cmat[k]);

    // Offsets into cmat (a triangular matrix)
    offset[k] = malloc(dim*sizeof(offset[k][0]));
    CHECKMEM(offset[k]);
    for(i=0; i < dim; ++i)
        offset[k][i] = ((2*dim - 1 - i)*i)/2;

    // Calculate matrices of eigenvectors in exact rational arithmetic
    Rational cvec[dim][dim], rvec[dim][dim], m;

    for(i=0; i<dim; ++i)
        for(j=0; j<dim; ++j)
            rvec[i][j] = cvec[i][j] = Rational_zero;
    
    for(j=2; j <= nLin; ++j) {
        jj = j-2;
        long lambda = -j*(j-1);
        cvec[jj][jj] = rvec[jj][jj] = Rational_unity;
        for(i=j-1; i > 1; --i) {
            ii = i-2;
            m = Rational_set(i*(i+1), i*(i-1) + lambda);
            cvec[ii][jj] = Rational_mul(cvec[ii+1][jj],m);
        }
        for(i=j+1; i <= nLin; ++i) {
            ii = i-2;
            m = Rational_set(i*(i-1), i*(i-1) + lambda);
            rvec[jj][ii] = Rational_mul(rvec[jj][ii-1], m);
        }
    }

    // Calculate coefficients of exponentials in x(t)
    // Convert to floating point and store in cmat[nLin]
    for(ii=0; ii<dim; ++ii) {
        for(jj=ii; jj<dim; ++jj) {
            cvec[ii][jj] = Rational_mul(cvec[ii][jj], rvec[jj][dim-1]);
            cmat[k][offset[k][ii] + jj] = Rational_ldbl(cvec[ii][jj]);
        }
    }

    // beta[i] is (i+2) choose 2
    Rational negBetaInv[dim];
    for(i=0; i<dim; ++i) {
        j = i+2;
        negBetaInv[i] = Rational_set(-2, j*(j-1));
    }

    Rational B[dim][dim+1];
    for(i=0; i<dim; ++i)
        B[i][0] = Rational_zero;
    B[dim-1][0] = Rational_set(-1,1);

    for(i=0; i < dim; ++i)
        for(j=0; j<dim; ++j)
            B[i][j+1] = cvec[i][j];

    for(i=dim-2; i>=0; --i)
        for(j=0; j<=dim; ++j)
            B[i][j] = Rational_add(B[i][j], B[i+1][j]);

    for(i=0; i<dim; ++i)
        for(j=0; j<=dim; ++j)
            B[i][j] = Rational_mul(B[i][j], negBetaInv[i]);

    bmat[k] = malloc(dim * (dim+1) * sizeof(**bmat));
    CHECKMEM(bmat[k]);

    // Convert B to floating point.
    for(i=0; i<dim; ++i)
        for(j=0; j<=dim; ++j)
            bmat[k][i*nLin + j] = Rational_ldbl(B[i][j]);
}

int MatCoal_initExterns(int n) {
    nsamples = n;
    int k;

    offset = malloc( (nsamples-2) * sizeof(offset[0]));
    CHECKMEM(offset);
        
    beta = malloc( (nsamples-2) * sizeof(beta[0]));
    CHECKMEM(beta);

    cmat = malloc( (nsamples-2) *sizeof(cmat[0]));
    CHECKMEM(cmat);

    for(k=2; k <= nsamples; ++k) {
        beta[k-2] = (k*(k-1))/2;
        init_dim(k-2);
    }

    return 0;
}

void MatCoal_freeExterns(void) {
    int i;
    if(nsamples == 0) {
        fprintf(stderr,"%s:%s:%d: can't free externs, because they"
                " aren't allocated.\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }
    
    for(i=2; i<nsamples; ++i) {
        free(cmat[i]);
        free(bmat[i]);
        free(offset[i]);
    }

    free(beta);
    free(cmat);
    free(bmat);
    free(offset);

    nsamples = 0;
}

/// Calculate the probability that, after v units of coalescent time,
/// there are 2,3,...(dim+1) lines of descent.
void MatCoal_project(int dim, double ans[dim], double v) {
    int i, j;
    double expn[dim];

    for(i=0; i<dim; ++i)
        expn[i] = exp(-v*beta[i]);

    // Multiply matrix cmat[dim-1] times vector expn
    for(i=0; i<dim; ++i) {
        ans[i] = 0.0;
        // Right-to-left sum accumulates small numbers first
        // to reduce error.
        for(j=dim-1; j >= i; --j)
            ans[i] += cmat[dim-1][offset[dim][i] + j] * expn[j];
    }
}

/// Vector of expected lengths of coalescent intervals during which
/// there were 2,3,...(dim+1) lines of descent. To get the expected
/// length of the interval with 1 line of descent, subtract the sum
/// of ans from v.
void MatCoal_ciLen(int dim, double ans[dim], double v) {
    int i, j;
    double expn[dim];

    for(i=0; i<dim; ++i)
        expn[i] = exp(-v*beta[i]);

    // Multiply matrix bmat[dim-1] times vector expn
    for(i=0; i<dim; ++i) {
        ans[i] = 0.0;
        for(j=dim-1; j >= i; --j)
            ans[i] += bmat[dim-1][i*(dim+1) + j] * expn[j];
        ans[i] += bmat[dim-1][i*(dim+1)];
    }
}
