#include "rational.h"
#include "matcoal.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

typedef struct MatCoal MatCoal;

static MatCoal * MatCoal_new(int nLin);
static void MatCoal_free(MatCoal *self);
static void MatCoal_print(MatCoal *self, FILE *fp);

// Data for an epoch with nLin lineages.
struct MatCoal {
    int nLin; // number of lineages in epoch

    // Array of dimension nLin-1. The i'th entry is (i+2) choose 2.
    // For example, if nLin=3, then beta has two entries: 2 choose 2
    // and 3 choose 2.  
    double *beta;
    
    /*
      Layout of upper triangular matrix gmat is row major:

      00 01 02 03 04  row 0: offset = 0
      xx 11 12 13 14  row 1: offset = dim-1
      xx xx 22 23 24  row 2: offset = 2*dim - 3
      xx xx xx 33 34  row 3: offset = 3*dim - 6
      xx xx xx xx 44  row 4: offset = 4*dim - 10

      offset[i] = i*dim - i*(i+1)/2 = ((2*dim - 1 - i)*i)/2
      x[i][j] is gmat[offset[i] + j]

      Number of stored elements is dim*(dim+1)/2. dim=nLin-1.
    */
    // Matrix of scaled column eigenvectors
    double *gmat;

    // Array of dim offsets, used to address elements of gmat.
    int *offset;

    // Matrix for calculating expected lengths of coalescent
    // intervals.  hmat is upper triangular and is stored just
    // like gmat. The ij-th element is hmat[offset[i] + j]
    double *hmat;

    // Vector for calculating expected lengths of coalescent intervals
    // Dimension is dim.
    double *z;
};

// Number of haploid samples in the data. 0 until initialized
static int nsamples = 0;

// Array of pointers to MatCoal objects. There are nsamples-1 entries
// in the array, and the i'th entry refers to an epoch that has
// i+2 lineages at the recent end of the epoch. The arrays for that
// epoch have i+1 rows.
static MatCoal **matcoal = NULL;

// Allocate and initialize an object of type MatCoal.
static MatCoal * MatCoal_new(int nLin) {
    long i, j, ii, jj, dim = nLin - 1;

    MatCoal *self = malloc(sizeof(MatCoal));
    CHECKMEM(self);
    memset(self, 0, sizeof(MatCoal));

    self->nLin = nLin;


    // i'th entry of beta is (i+2) choose 2
    self->beta = malloc(dim * sizeof(self->beta[0]));
    for(i=0; i<dim; ++i)
        self->beta[i] = ((i+2)*(i+1))/2;

    size_t size = sizeof(self->gmat[0]) * (dim*(dim+1))/2;
    self->gmat = malloc(size);
    CHECKMEM(self->gmat);

    // Offsets into gmat (a triangular matrix)
    self->offset = malloc( dim * sizeof(self->offset[0]));
    CHECKMEM(self->offset);
    for(i=0; i < dim; ++i)
        self->offset[i] = ((2*dim - 1 - i)*i)/2;
    
    // Calculate matrices of eigenvectors in exact rational arithmetic
    Rational cvec[dim][dim], rvec[dim][dim], m;

    for(i=0; i<dim; ++i)
        for(j=0; j<dim; ++j)
            rvec[i][j] = cvec[i][j] = Rational_zero;
    
    for(j=2; j <= self->nLin; ++j) {
        jj = j-2;
        long lambda = -j*(j-1);
        cvec[jj][jj] = rvec[jj][jj] = Rational_unity;
        for(i=j-1; i > 1; --i) {
            ii = i-2;
            m = Rational_set(i*(i+1), i*(i-1) + lambda);
            cvec[ii][jj] = Rational_mul(cvec[ii+1][jj],m);
        }
        for(i=j+1; i <= self->nLin; ++i) {
            ii = i-2;
            m = Rational_set(i*(i-1), i*(i-1) + lambda);
            rvec[jj][ii] = Rational_mul(rvec[jj][ii-1], m);
        }
    }

    // Calculate coefficients of exponentials in x(t)
    // Convert to floating point and store in gmat.
    for(ii=0; ii<dim; ++ii) {
        for(jj=ii; jj<dim; ++jj) {
            cvec[ii][jj] = Rational_mul(cvec[ii][jj], rvec[jj][dim-1]);
            self->gmat[self->offset[ii] + jj] = Rational_ldbl(cvec[ii][jj]);
        }
    }
    
    // beta[i] is (i+2) choose 2
    // negBetaInv[i] is -1/beta[i]
    Rational negBetaInv[dim];
    for(i=0; i<dim; ++i) {
        j = i+2;
        negBetaInv[i] = Rational_set(-2, j*(j-1));
    }

    Rational H[dim][dim];

    // Initial H is cvec
    for(i=0; i < dim; ++i)
        for(j=0; j<dim; ++j)
            H[i][j] = cvec[i][j];

    // Cumulative sum, so that i'th row is sum of i..(dim-1) initial
    // rows.
    for(i=dim-2; i>=0; --i)
        for(j=0; j<dim; ++j)
            H[i][j] = Rational_add(H[i][j], H[i+1][j]);

    // Weight row i by negBetaInv[i]
    for(i=0; i<dim; ++i)
        for(j=0; j<dim; ++j)
            H[i][j] = Rational_mul(H[i][j], negBetaInv[i]);

    self->hmat = malloc(size);
    CHECKMEM(self->hmat);
    self->z = malloc(dim * sizeof(self->z[0]));
    CHECKMEM(self->z);

    // Convert to floating point and store in z and hmat.
    for(ii=0; ii<dim; ++ii) {
        self->z[ii] = Rational_ldbl(negBetaInv[ii]);
        for(jj=ii; jj<dim; ++jj)
            self->hmat[self->offset[ii] + jj] = Rational_ldbl(H[ii][jj]);
    }
    return self;
}

void MatCoal_initExterns(long nsamp) {
    matcoal = malloc((nsamp-1) * sizeof(matcoal[0]));
    CHECKMEM(matcoal);

    nsamples = nsamp;

    for(long i=2; i <= nsamp; ++i)
        matcoal[i-2] = MatCoal_new(i);
}

static void MatCoal_free(MatCoal *self) {
    free(self->beta);
    free(self->gmat);
    free(self->offset);
    free(self->z);
    free(self->hmat);
    free(self);
}

void MatCoal_freeExterns(void) {
    if(nsamples == 0) {
        fprintf(stderr,"%s:%s:%d: can't free externs, because they"
                " aren't allocated.\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }
    
    for(int i=0; i < nsamples-1; ++i)
        MatCoal_free(matcoal[i]);

    free(matcoal);

    nsamples = 0;
}

/// Calculate the probability that, after v units of coalescent time,
/// there are 2,3,...(dim+1) lines of descent.
void MatCoal_project(int dim, double ans[dim], double v) {
    int i, j;
    int ndx = dim-1;
    MatCoal *mc = matcoal[ndx];
    assert(dim == mc->nLin - 1);
    
    double expn[dim];

    for(i=0; i<dim; ++i)
        expn[i] = exp(-v*mc->beta[i]);

    // Multiply matrix gmat[dim-1] times vector expn
    for(i=0; i<dim; ++i) {
        ans[i] = 0.0;
        // Right-to-left sum accumulates small numbers first
        // to reduce error.
        for(j=dim-1; j >= i; --j)
            ans[i] += mc->gmat[mc->offset[i] + j] * expn[j];
    }
}

/// Vector of expected lengths of coalescent intervals during which
/// there were 2,3,...(dim+1) lines of descent. To get the expected
/// length of the interval with 1 line of descent, subtract the sum
/// of ans from v.
void MatCoal_ciLen(int dim, double ans[dim], double v) {
    int i, j;
    int ndx = dim-1;
    MatCoal *mc = matcoal[ndx];
    assert(dim == mc->nLin - 1);

    double expn[dim];

    for(i=0; i<dim; ++i)
        expn[i] = exp(-v*mc->beta[i]);

    // ans = z + H*expn
    for(i=0; i<dim; ++i) {
        ans[i] = 0.0;
        for(j=dim-1; j >= i; --j)
            ans[i] += mc->hmat[mc->offset[i] + j] * expn[j];
        ans[i] += mc->z[i];
    }
}

void MatCoal_printAll(FILE *fp) {
    fprintf(fp, "nsamples=%d\n", nsamples);
    for(int i=0; i < nsamples-1; ++i)
        MatCoal_print(matcoal[i], fp);
}


// Print print a MatCoal object
static void MatCoal_print(MatCoal *self, FILE *fp) {
    int i, j, dim = self->nLin-1;

    fprintf(fp, "\nnLin=%d\n", self->nLin);
    fprintf(fp, "beta:");
    for(i=0; i < dim; ++i)
        fprintf(fp, " %lf", self->beta[i]);
    putc('\n', fp);

    fprintf(fp, "G:\n");
    for(i=0; i<dim; ++i) {
        for(j=0; j<i; ++j)
            fprintf(fp," 0");
        for(j=i; j<dim; ++j)
            fprintf(fp," %lf", self->gmat[self->offset[i] + j]);
        putc('\n', fp);
    }

    fprintf(fp, "z:");
    for(i=0; i<dim; ++i)
        fprintf(fp, " %lf", self->z[i]);
    putc('\n', fp);

    fprintf(fp, "H:\n");
    for(i=0; i<dim; ++i) {
        for(j=0; j<i; ++j)
            fprintf(fp, " 0");
        for(j=i; j<dim; ++j)
            fprintf(fp," %lf", self->hmat[self->offset[i] + j]);
        putc('\n', fp);
    }
}
