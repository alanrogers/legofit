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
// cmat[k] is a matrix of dimension kXk. The first two entries,
// cmat[0] and cmat[1], are NULL.
static double **cmat = NULL;

// offset[k] is a pointer to an array of k offsets, which are
// used to address elements of the triangular matrix cmat[k]. The
// first two entries of offset, which correspond to nonexistent
// matrices of dimension 0 and 1, are NULL
static int **offset;

// Array of pointers to matrices for calculating expected branch
// lengths. First two pointers are NULL. For k>1, kth matrix is
// rectangular, with dimension (k-1) X k, but allocated as a linear
// array. (i,j)th element of k'th matrix is bmat[k][i*k + j]
static double **bmat = NULL;

// beta[i] = (i+1)*(i+2)/2, i=0..(nsamples-2)
static double *beta = NULL;

// Initialize arrays for epochs in which there are nLin Lineages.
static void init_dim(int nLin) {
    int i, j, ii, jj;
    const int dim = nLin-1;
    size_t size = sizeof(cmat[nLin][0]) * (dim*(dim+1))/2;
    cmat[nLin] = malloc(size);
    CHECKMEM(cmat[nLin]);

    // Offsets into cmat (a triangular matrix)
    offset[nLin] = malloc(dim*sizeof(offset[nLin][0]));
    CHECKMEM(offset[nLin]);
    for(i=0; i < dim; ++i)
        offset[nLin][i] = ((2*dim - 1 - i)*i)/2;

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
            cmat[nLin][offset[nLin][ii] + jj] = Rational_ldbl(cvec[ii][jj]);
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

    bmat[nLin] = malloc(dim * (dim+1) * sizeof(**bmat));
    CHECKMEM(bmat[nLin]);

    // Convert B to floating point.
    for(i=0; i<dim; ++i)
        for(j=0; j<=dim; ++j)
            bmat[nLin][i*nLin + j] = Rational_ldbl(B[i][j]);
}

int MatCoal_initExterns(int n) {
    nsamples = n;
    int k;

    offset = malloc( (nsamples+1) * sizeof(offset[0]));
    CHECKMEM(offset);
        
    beta = malloc( (nsamples-1) * sizeof(beta[0]));
    CHECKMEM(beta);

    cmat = malloc( (nsamples+1) *sizeof(cmat[0]));
    CHECKMEM(cmat);

    cmat[0] = cmat[1] = NULL;
    offset[0] = offset[1] = NULL;
    for(k=2; k <= nsamples; ++k) {
        beta[k-2] = (k*(k-1))/2;

        init_dim(k);
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


int MatCoal_project(int dim, double ans[dim], double v) {
    int i, j;
    double expn[dim];

    for(i=0; i<dim; ++i)
        expn[i] = exp(-v*beta[i]);

    // Multiply matrix cmat[dim+2] times vector expn
    for(i=0; i<dim; ++i) {
        ans[i] = 0.0;
        // right-to-left sum accumulates small numbers first
        for(j=dim-1; j >= i; --j)
            ans[i] += cmat[dim+2][offset[dim+2] + j] * expn[j];
    }
}
