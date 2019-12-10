/**
 * @file mpfrmatcoal.c
 * @brief Spectral decomposition of matrix coalescent in extended precision.
 * @internal
 * Copyright (C) 2015,2019 Alan R. Rogers
 *
 * @copyright Copyright (c) 2015, 2019, Alan R. Rogers
 * This file is released under the Internet Systems Consortium
 * License, which can be found in file "LICENSE".
 *
 * Alan R. Rogers, Department of Anthropology, University of Utah,
 * Salt Lake City, UT 84112. Email: rogers at anthro.utah.edu
 **/

#include "mpfrmatcoal.h"
#include "misc.h"
#include <string.h>
#include <stdarg.h>
#include <mpfr.h>
#include <math.h>

const mpfr_rnd_t rnd = MPFR_RNDN;  // round to nearest
const mpfr_prec_t precision = 256; // bits of precision

/*
   0  1  2  3  4  col j
   0  1  3  6 10  offset = j*(j+1)/2

  00 01 02 03 04  
  xx 11 12 13 14  
  xx xx 22 23 24  
  xx xx xx 33 34  
  xx xx xx xx 44  

  offset[j] = j*(j+1)/2
  x[i][j] is array[offset[j] + i]
 */
struct MpfrMatCoal {
    unsigned nSamples;
    unsigned dim; // dimension of matrices and vectors
    unsigned nPairs;
    mpfr_t *rvec; // row eigenvectors
    mpfr_t *cvec; // column eigenvectors
    mpfr_t *beta; // beta[i-2] = i*(i-1)/2
    mpfr_t *lambda; // eigenvalues
    mpfr_t *work, *work2, v, w, x, y, z, dv; // temporaries
    unsigned *offset;  // offset[j] = j*(j+1)/2
};

void prUTMat(unsigned dim, mpfr_t *mat, unsigned prWid, unsigned offset[dim]);
static void print_mpfr_vec(unsigned n, mpfr_t v[n]);

/// Allocate and initialize a new vector of mpfr_t values.
/// If x==NULL on input, the vector is allocated but not initialized.
MpfrVec *MpfrVec_new(unsigned dim, long double x[dim]) {
    MpfrVec *new = malloc(sizeof(new[0]));
    CHECKMEM(new);

    new->dim = dim;
    new->x = malloc(dim * sizeof(new->x[0]));
    CHECKMEM(new->x);
    unsigned i;
    for(i=0; i < dim; ++i)
        mpfr_init2(new->x[i], precision);

    if(x != NULL)
        for(i=0; i<dim; ++i)
            mpfr_set_ld(new->x[i], x[i], rnd);

    return new;
}

/// Deallocate MpfrVec
void MpfrVec_free(MpfrVec *self) {
    unsigned i;
    for(i=0; i < self->dim; ++i)
        mpfr_clear(self->x[i]);
    free(self->x);
    free(self);
}

/// Copy contents of MpfrVec into array x.
void MpfrVec_get(MpfrVec *self, unsigned dim, long double x[dim]) {
    assert(dim == self->dim);
    unsigned i;

    for(i=0; i < dim; ++i)
        x[i] = mpfr_get_ld(self->x[i], rnd);
}

/// Set MpfrVec equal to array x.
void MpfrVec_set(MpfrVec *self, unsigned dim, const long double x[dim]) {
    assert(dim == self->dim);
    unsigned i;

    for(i=0; i < dim; ++i)
        mpfr_set_ld(self->x[i], x[i], rnd);
}

void MpfrVec_print(MpfrVec *self, FILE *ofp) {
    unsigned i;
    for(i=0; i < self->dim; ++i) {
        fprintf(ofp,"%3u ", i);
        mpfr_fprintf(ofp, "%Rf\n", self->x[i], rnd);
    }
}

static void print_mpfr_vec(unsigned n, mpfr_t v[n]) {
    unsigned i;

    for(i=0; i<n; ++i)
        mpfr_printf(" %Rf", v[i]);
    putchar('\n');
    
}

MpfrMatCoal *MpfrMatCoal_new(unsigned nSamples) {
    long i, j, ii, jj;
    MpfrMatCoal *self = malloc( sizeof *self );
    CHECKMEM(self);

    self->nSamples = nSamples;
    self->dim = nSamples-1;
    self->nPairs = (self->dim * (self->dim+1))/2;

    mpfr_inits2(precision, self->v, self->w, self->x, self->y,
                self->z, self->dv, (mpfr_ptr) 0);

    self->offset = malloc(self->dim * sizeof self->offset[0]);
    CHECKMEM(self->offset);

    self->rvec = malloc(self->nPairs * sizeof self->rvec[0]);
    CHECKMEM(self->rvec);

    self->cvec = malloc(self->nPairs * sizeof self->cvec[0]);
    CHECKMEM(self->cvec);

    for(i=0; i < self->nPairs; ++i) {
        mpfr_init2(self->rvec[i], precision);
        mpfr_init2(self->cvec[i], precision);
    }

    for(i=0; i < self->dim; ++i)
        self->offset[i] = (i*(i + 1))/2;

    for(j=2; j <= nSamples; ++j) {
        jj = j-2;
        mpfr_set_si(self->cvec[jj + self->offset[jj]], 1L, rnd);
        mpfr_set_si(self->rvec[jj + self->offset[jj]], 1L, rnd);
        for(i = j-1; i > 1; --i) {
            ii = i-2;
            mpfr_set_si(self->x, i*(i+1L), rnd);
            mpfr_set_si(self->y, i*(i-1L) - j*(j-1L), rnd);
            mpfr_div(self->z, self->x, self->y, rnd);
            // now z = i*(i+1)/(i*(i-1) - j*(j-1))
            mpfr_mul(self->cvec[ii + self->offset[jj]],
                     self->cvec[ii+1 + self->offset[jj]],
                     self->z, rnd);
        }
        for(i=j+1; i <= nSamples; ++i) {
            ii = i-2;
            mpfr_set_si(self->x, i*(i-1L), rnd);
            mpfr_set_si(self->y, i*(i-1L) - j*(j-1L), rnd);
            mpfr_div(self->z, self->x, self->y, rnd);
            // z = i*(i-1)/(i*(i-1) - j*(j-1))
            mpfr_mul(self->rvec[jj+self->offset[ii]],
                     self->rvec[jj+self->offset[ii-1]], self->z, rnd);
        }
    }

    self->beta = malloc(self->dim * sizeof self->beta[0]);
    CHECKMEM(self->beta);

    self->lambda = malloc(self->dim * sizeof self->lambda[0]);
    CHECKMEM(self->lambda);

    self->work = malloc(self->dim * sizeof self->work[0]);
    CHECKMEM(self->work);

    self->work2 = malloc(self->dim * sizeof self->work2[0]);
    CHECKMEM(self->work2);

    for(i=0; i < self->dim; ++i) {
        j = i+2;
        mpfr_init2(self->beta[i], precision);
        mpfr_init2(self->lambda[i], precision);
        mpfr_init2(self->work[i], precision);
        mpfr_init2(self->work2[i], precision);
        mpfr_set_si(self->beta[i], (j*(j-1L))/2L, rnd);
    }
    return self;
}

void MpfrMatCoal_free(MpfrMatCoal *self) {
    unsigned i;

    mpfr_clears(self->v, self->w, self->x, self->y, self->z, self->dv,
                (mpfr_ptr) 0);

    for(i=0; i < self->nPairs; ++i)
        mpfr_clears(self->rvec[i], self->cvec[i], (mpfr_ptr) 0);

    for(i=0; i < self->dim; ++i)
        mpfr_clears(self->beta[i], self->lambda[i], self->work[i],
                    self->work2[i], (mpfr_ptr) 0);

    free(self->beta);
    free(self->cvec);
    free(self->rvec);
    free(self->offset);
    free(self);
}

// Project vector an initial vector, (0,0,...,0,1)', forward by v time
// units. Put result into x.
void MpfrMatCoal_project(MpfrMatCoal *self, int dim, long double ans[dim],
                         long double v) {
    if(dim != self->dim) {
        fprintf(stderr,"%s:%d: dimension mismatch\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    int i;
    
    for(i=0; i < self->dim; ++i){
        // Set work equal to right-most column of rvec. This is the
        // same as left-multiplying initial vector by rvec.
        mpfr_set(self->work[i], self->rvec[i + self->offset[self->dim-1]],
                 rnd);
    }

    // Left-multiply resulting vector by diag(exp(-beta[i]*v))
    for(i=0; i < self->dim; ++i) {
        // self->x = exp(-beta[i]*v)
        mpfr_set_ld(self->x, v, rnd);
        mpfr_mul(self->x, self->x, self->beta[i], rnd); 
        mpfr_neg(self->x, self->x, rnd); 
        mpfr_exp(self->x, self->x, rnd); 
        mpfr_mul(self->work[i], self->work[i], self->x, rnd);
    }

    // Left-multiply by matrix of column eigenvectors.
    UTmatXvec(self->dim, self->work2, self->cvec, self->offset, self->work);

    for(i=0; i < dim; ++i)
        ans[i] = mpfr_get_ld(self->work2[i], rnd);
    
}

/// Form matrix product y = A*x, where y is a vector, A an
/// upper-triangular matrix, and x is a vector.  x and y have
/// dimension dim, and A has dimension dim X dim. offset[j] is the
/// offset of the beginning of the j'th column within matrix A.
/// A should be laid out so that element (i,j), i.e. row i and column
/// j, is at position A[i + offset[j]].
void UTmatXvec(unsigned dim, mpfr_t *y, mpfr_t *A, unsigned *offset,
               mpfr_t *x) {
    unsigned i, j;
    mpfr_t tmp;
    mpfr_init2(tmp, precision);

    for(i=0; i<dim; ++i) {
        mpfr_set_si(y[i], 0, rnd);
        for(j=i; j<dim; ++j) {
            mpfr_mul(tmp, A[i + offset[j]], x[j], rnd);
            mpfr_add(y[i], y[i], tmp, rnd);
        }
    }

    mpfr_clear(tmp);
}

void MpfrMatCoal_print(MpfrMatCoal *self) {
    printf("nSamples=%u nPairs=%d precision=%ld\n",
           self->nSamples, self->nPairs, precision);

    long i;

    printf("beta:");
    for(i=0; i < self->dim; ++i)
        mpfr_printf(" %RNf", self->beta[i]);
    putchar('\n');
    
    printf("Row eigenvectors:\n");
    prUTMat(self->dim, self->rvec, 17, self->offset);

    printf("Column eigenvectors:\n");
    prUTMat(self->dim, self->cvec, 17, self->offset);
}

void MpfrMatCoal_ciLen(MpfrMatCoal *self, unsigned dim, long double m[dim],
                       long double v) {
    assert(self);

    if(dim != self->dim)
        eprintf("%s:%s:%d: dimension mismatch. dim=%u but self->dim=%u\n",
                __FILE__,__func__,__LINE__, dim, self->dim);

    unsigned i, j;

    for(j=0; j<dim; ++j)
        mpfr_set_si(self->lambda[j], 0, rnd);

    // v = t/twoN
    mpfr_set_ld(self->v, v, rnd);

    for(j=0; j<dim; ++j) {
        // z: expm1(-beta*v)
        mpfr_mul(self->z, self->beta[j], self->v, rnd); // beta*v
        mpfr_neg(self->z, self->z, rnd);           // -beta*v
        mpfr_expm1(self->z, self->z, rnd);         // expm1(-beta*v)

        // lambda[j] += expm1(-beta*v)
        mpfr_add(self->lambda[j], self->lambda[j], self->z, rnd);
    }

    // lambda[i] /= -beta[i]. Equivalent to multiplying by inv(B).
    for(i=0; i < dim; ++i){
        mpfr_div(self->lambda[i], self->lambda[i], self->beta[i], rnd);
        mpfr_neg(self->lambda[i], self->lambda[i], rnd);
    }

    // Matrix multiplication
    for(i=0; i < dim; ++i){
        // Set work equal to right-most column of rvec. This is the same
        // as right-multiplying rvec by (0,0,...,1)'.
        mpfr_set(self->work[i], self->rvec[i + self->offset[dim-1]], rnd);

        // Multiply elements by lambda: same as left-multiplying
        // work by diag(lambda).  
        mpfr_mul(self->work[i], self->work[i], self->lambda[i], rnd);
    }

    // left-multiply work by cvec; result in work2
    UTmatXvec(dim, self->work2, self->cvec, self->offset, self->work);

    for(i=0; i < dim; ++i)
        m[i] = mpfr_get_ld(self->work2[i], rnd);
}

void prUTMat(unsigned dim, mpfr_t *mat, unsigned prWid, unsigned offset[dim]) {
    int i, j;
    for(i=0; i<dim; ++i) {
        for(j=0; j < i; ++j) 
            printf(" %*d", prWid, 0);
        for(j=i; j < dim; ++j) {
            putchar(' ');
#if 0
            cWritten = mpfr_out_str(stdout, 10, prWid,
                                    mat[i + offset[j]],
                                    rnd);
#else
            mpfr_printf("%*.6RNf", prWid,mat[i + offset[j]]);
#endif
        }
        putchar('\n');
    }
}

#ifdef TEST

#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {

    int verbose=0, ok=1;
    long double err, maxerr, errTol = DBL_EPSILON;

    if(argc > 1) {
        if(argc!=2 || 0!=strcmp(argv[1], "-v")) {
            fprintf(stderr,"usage: xmpfrmatcoal [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    unsigned i, nSamples = 10;
    MpfrMatCoal *mcs = MpfrMatCoal_new(nSamples);

    if(verbose)
        MpfrMatCoal_print(mcs);

    unsigned  dim = nSamples-1;
    long double x[dim], y[dim], m[dim];
    memset(x, 0, sizeof(x));
    x[dim-1] = 1.0L;

    // project 0 units into future */
    long double v = 0.0L;
    MpfrMatCoal_project(mcs, dim, y, v);
    maxerr = 0.0L;
    if(verbose)
        printf("v=%Lf\n", v);
    for(i=0; i < dim; ++i) {
        err = fabsl(y[i] - x[i]);  // y should equal x
        maxerr = fmaxl(maxerr, err);
        if(verbose)
            printf("%3d %10.8Lf -> %12.8Lf\n", i+2, x[i], y[i]);
    }

    if(verbose) {
        printf("After _project(0):\n");
        for(i=0; i<dim; ++i)
            printf("y[%d] = %Lf\n", i, y[i]);
    }

    if(maxerr > errTol) {
        ok = 0;
        fprintf(stderr,"%s:%s:%d: maxerr=%Le FAIL\n",
                __FILE__, __func__, __LINE__, maxerr);
    }

    // project 1 unit forward
    v = 1.0L;
    MpfrMatCoal_project(mcs, dim, y, v);

    if(verbose) {
        printf("After _project(1):\n");
        for(i=0; i<dim; ++i)
            printf("y[%d] = %Lf\n", i, y[i]);
    }

    // When v is very large, vec should be zero.
    v = INFINITY;
    maxerr = 0.0L;
    MpfrMatCoal_project(mcs, dim, y, v);
    for(i=0; i < dim; ++i) {
        err = fabsl(y[i]);  // should equal 0
        maxerr = fmaxl(maxerr, err);
        if(verbose)
            printf("%3d %10.8Lf -> %12.8Lf\n", i+2, 
                   (i==dim-1 ? 1.0L : 0.0L), y[i]);
    }

    if(verbose) {
        printf("After _project(inf):\n");
        for(i=0; i<dim; ++i)
            printf("y[%d] = %Lf\n", i, y[i]);
    }

    if(maxerr > errTol) {
        ok = 0;
        fprintf(stderr,"%s:%s:%d: maxerr=%Le FAIL\n",
                __FILE__, __func__, __LINE__, maxerr);
    }

    unitTstResult("MpfrMatCoal_project", ok ? "OK" : "FAIL");

    ok = 1;

    v = 0.0L;

    MpfrMatCoal_ciLen(mcs, dim, m, v);

    if(verbose) {
        printf("After _ciLen(0):\n");
        for(i=0; i<dim; ++i)
            printf("m[%d] = %Lf\n", i, m[i]);
    }

    unitTstResult("MpfrMatCoal_ciLen", ok ? "OK" : "FAIL");
    return 0;
}
#endif
