#include "rational.h"
#include <stdio.h>
#include <string.h>

// Number of haploid samples in the data. 0 until initialized
static int nsamples = 0;

/*
  Layout of upper triangular matrices is row major:

  00 01 02 03 04  row 0: offset = 0
  xx 11 12 13 14  row 1: offset = dim
  xx xx 22 23 24  row 2: offset = 2*dim - 1
  xx xx xx 33 34  row 3: offset = 3*dim - 3
  xx xx xx xx 44  row 4: offset = 4*dim - 6

  offset[i] = i*dim - i*(i-1)/2 = ((2*dim+1 - i)*i)/2
  x[i][j] is array[offset[i] + j - i]

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

// Array of pointers to matrices for calculating expected
// branch lengths. Each matrix is rectangular.
static double **bmat = NULL;

// beta[i] = i*(i-1)/2. beta[0] and beta[1] are not used.
static double *beta = NULL;

int MatCoal_initExterns(int n) {
    nsamples = n;
    int i, j, k;

    offset = malloc( (nsamples+1) * sizeof(offset[0]));
    CHECKMEM(offset);
        
    beta = malloc( (nsamples+1) * sizeof(beta[0]));
    CHECKMEM(beta);

    cmat = malloc( (nsamples+1) *sizeof(cvec[0]));
    CHECKMEM(cmat);

    cmat[0] = cmat[1] = NULL;
    offset[0] = offset[1] = NULL;
    beta[0]=beta[1]=0.0;
    for(k=2; k <= nsamples; ++k) {
        beta[k] = (k*(k-1))/2;

        offset[k] = malloc(k*sizeof(offset[k][0]));
        CHECKMEM(offset[k]);
        for(i=0; i < k; ++i)
            offset[k][i] = ((2*k+1 - i)*i)/2;

        cmat[k] = new_cmat(k);
    }

    return 0;
}
