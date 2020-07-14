#include "rational.h"
#include "mpfrmatcoal.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpfr.h>

typedef struct MpfrMatCoal MpfrMatCoal;

static MpfrMatCoal * MpfrMatCoal_new(int nLin);
static void MpfrMatCoal_free(MpfrMatCoal *self);
static void MpfrMatCoal_print(MpfrMatCoal *self, FILE *fp);

const mpfr_rnd_t rnd = MPFR_RNDN;  // round to nearest
const mpfr_prec_t precision = 256; // bits of precision

// Data for an epoch with nLin lineages.
struct MpfrMatCoal {
    int nLin; // number of lineages in epoch

    // Array of dimension nLin-1. The i'th entry is (i+2) choose 2.
    // For example, if nLin=3, then beta has two entries: 2 choose 2
    // and 3 choose 2.  
    mpfr_t *beta;
    
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
    mpfr_t *gmat;

    // Array of dim offsets, used to address elements of gmat.
    int *offset;

    // Matrix for calculating expected lengths of coalescent
    // intervals.  hmat is upper triangular and is stored just
    // like gmat. The ij-th element is hmat[offset[i] + j]
    mpfr_t *hmat;
    
    // Vector for calculating expected lengths of coalescent intervals
    // Dimension is dim.
    mpfr_t *z;

    mpfr_t *expn; // temporary storage
    mpfr_t v, x, y;
};

// Number of haploid samples in the data. 0 until initialized
static int nsamples = 0;

// Array of pointers to MpfrMatCoal objects. There are nsamples-1 entries
// in the array, and the i'th entry refers to an epoch that has
// i+2 lineages at the recent end of the epoch. The arrays for that
// epoch have i+1 rows.
static MpfrMatCoal **matcoal = NULL;

// Allocate and initialize an object of type MpfrMatCoal.
static MpfrMatCoal * MpfrMatCoal_new(int nLin) {
    long i, j, ii, jj, dim = nLin - 1;

    MpfrMatCoal *self = malloc(sizeof(MpfrMatCoal));
    CHECKMEM(self);
    memset(self, 0, sizeof(MpfrMatCoal));

    self->nLin = nLin;

    mpfr_init2(self->v, precision);
    mpfr_init2(self->x, precision);
    mpfr_init2(self->y, precision);

    // i'th entry of beta is (i+2) choose 2
    self->beta = malloc(dim * sizeof(self->beta[0]));
    for(i=0; i<dim; ++i) {
        mpfr_init2(self->beta[i], precision);
        mpfr_set_si(self->beta[i], ((i+2)*(i+1))/2, rnd);
    }

    // temporary storage
    self->expn = malloc(dim * sizeof(self->expn[0]));
    for(i=0; i<dim; ++i)
        mpfr_init2(self->expn[i], precision);
    
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

    int npairs = (dim*(dim+1))/2;
    size_t size = sizeof(self->gmat[0]) * npairs;
    self->gmat = malloc(size);
    CHECKMEM(self->gmat);

    for(i=0; i<npairs; ++i)
        mpfr_init2(self->gmat[i], precision);

    mpfr_t num, den;
    mpfr_init2(num, precision);
    mpfr_init2(den, precision);

    // Offsets into gmat (a triangular matrix)
    self->offset = malloc( dim * sizeof(self->offset[0]));
    CHECKMEM(self->offset);
    for(i=0; i < dim; ++i)
        self->offset[i] = ((2*dim - 1 - i)*i)/2;

    // Calculate coefficients of exponentials in x(t)
    // Convert to floating point and store in gmat.
    for(ii=0; ii<dim; ++ii) {
        for(jj=ii; jj<dim; ++jj) {
            cvec[ii][jj] = Rational_mul(cvec[ii][jj], rvec[jj][dim-1]);
            mpfr_set_si(num, cvec[ii][jj].num, rnd); // numerator
            mpfr_set_si(den, cvec[ii][jj].den, rnd); // denominator
            mpfr_div(self->gmat[self->offset[ii] + jj], num, den, rnd); // ratio
        }
    }
    
    // beta[i] is (i+2) choose 2
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

    for(i=0; i < npairs; ++i)
        mpfr_init2(self->hmat[i], precision);

    for(i=0; i < dim; ++i)
        mpfr_init2(self->z[i], precision);

    // Convert to floating point and store in z and hmat.
    for(i=0; i<dim; ++i) {
        mpfr_set_si(num, negBetaInv[i].num, rnd); // numerator
        mpfr_set_si(den, negBetaInv[i].den, rnd); // denominator
        mpfr_div(self->z[i], num, den, rnd); // ratio
        mpfr_neg(self->z[i], self->z[i], rnd); // negate
        for(j=i; j<dim; ++j) {
            mpfr_set_si(num, H[i][j].num, rnd); // numerator
            mpfr_set_si(den, H[i][j].den, rnd); // denominator
            mpfr_div(self->hmat[self->offset[i] + j], num, den, rnd); // ratio
        }
    }

    mpfr_clear(num);
    mpfr_clear(den);

    return self;
}

void MpfrMatCoal_initExterns(long nsamp) {
    matcoal = malloc((nsamp-1) * sizeof(matcoal[0]));
    CHECKMEM(matcoal);

    nsamples = nsamp;

    for(long i=2; i <= nsamp; ++i)
        matcoal[i-2] = MpfrMatCoal_new(i);
}

static void MpfrMatCoal_free(MpfrMatCoal *self) {

    int i, dim = self->nLin - 1;

    mpfr_clear(self->v);
    mpfr_clear(self->x);
    mpfr_clear(self->y);
    
    for(i=0; i<dim; ++i) {
        mpfr_clear(self->beta[i]);
        mpfr_clear(self->z[i]);
    }

    int npairs = (dim*(dim+1))/2;

    for(i=0; i < npairs; ++i)
        mpfr_clear(self->gmat[i]);

    for(i=0; i < npairs; ++i)
        mpfr_clear(self->hmat[i]);

    free(self->offset);
    free(self->beta);
    free(self->gmat);
    free(self->hmat);
    free(self->z);
    free(self);
}

void MpfrMatCoal_freeExterns(void) {
    if(nsamples == 0) {
        fprintf(stderr,"%s:%s:%d: can't free externs, because they"
                " aren't allocated.\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }
    
    for(int i=0; i < nsamples-1; ++i)
        MpfrMatCoal_free(matcoal[i]);

    free(matcoal);

    nsamples = 0;
}

/// Calculate the probability that, after v units of coalescent time,
/// there are 2,3,...(dim+1) lines of descent.
void MpfrMatCoal_project(int dim, double ans[dim], double v) {
    int i, j;
    int ndx = dim-1;
    MpfrMatCoal *mc = matcoal[ndx];
    assert(dim == mc->nLin - 1);

    mpfr_set_d(mc->v, v, rnd);

    // set expn[i] = exp(-beta[i]*v)
    for(i=0; i<dim; ++i) {
        mpfr_mul(mc->x, mc->v, mc->beta[i], rnd); // x = beta[i]*v
        mpfr_neg(mc->x, mc->x, rnd);              // x = -beta[i]*v
        mpfr_exp(mc->expn[i], mc->x, rnd);        // expn[i] = exp(-beta[i]*v)
    }

    // Multiply matrix gmat[dim-1] times vector expn
    for(i=0; i<dim; ++i) {
        mpfr_set_d(mc->x, 0.0, rnd);                       // x = 0
        // Right-to-left sum accumulates small numbers first
        // to reduce error.
        for(j=dim-1; j >= i; --j) {
            mpfr_mul(mc->y, mc->gmat[mc->offset[i] + j], mc->expn[j], rnd);
            mpfr_add(mc->x, mc->x, mc->y, rnd);
        }
        ans[i] = mpfr_get_d(mc->x, rnd);
    }
}

/// Vector of expected lengths of coalescent intervals during which
/// there were 2,3,...(dim+1) lines of descent. To get the expected
/// length of the interval with 1 line of descent, subtract the sum
/// of ans from v.
void MpfrMatCoal_ciLen(int dim, double ans[dim], double v) {
    int i, j;
    int ndx = dim-1;
    MpfrMatCoal *mc = matcoal[ndx];
    assert(dim == mc->nLin - 1);

    mpfr_set_d(mc->v, v, rnd);

    // set expn[i] = exp(-beta[i]*v)
    for(i=0; i<dim; ++i) {
        mpfr_mul(mc->x, mc->v, mc->beta[i], rnd); // x = beta[i]*v
        mpfr_neg(mc->x, mc->x, rnd);              // x = -beta[i]*v
        mpfr_exp(mc->expn[i], mc->x, rnd);        // expn[i] = exp(-beta[i]*v)
    }

    // ans = z + H*expn
    for(i=0; i<dim; ++i) {
        mpfr_set_d(mc->x, 0.0, rnd);                       // x = 0
        for(j=dim-1; j >= i; --j) {
            mpfr_mul(mc->y, mc->hmat[mc->offset[i] + j], mc->expn[j], rnd);
            mpfr_add(mc->x, mc->x, mc->y, rnd);
        }
        mpfr_add(mc->x, mc->x, mc->z[i], rnd);
        ans[i] = mpfr_get_d(mc->x, rnd);
    }
}

void MpfrMatCoal_printAll(FILE *fp) {
    fprintf(fp, "nsamples=%d\n", nsamples);
    for(int i=0; i < nsamples-1; ++i)
        MpfrMatCoal_print(matcoal[i], fp);
}


// Print print a MpfrMatCoal object
static void MpfrMatCoal_print(MpfrMatCoal *self, FILE *fp) {
    int i, j, dim = self->nLin-1;

    fprintf(fp, "\nMPFR: nLin=%d\n", self->nLin);
    fprintf(fp, "MPFR: beta:");
    for(i=0; i < dim; ++i)
        mpfr_fprintf(fp, " %Rf", self->beta[i], rnd);
    putc('\n', fp);

    fprintf(fp, "MPFR: G:\n");
    for(i=0; i<dim; ++i) {
        for(j=0; j<i; ++j)
            fprintf(fp," 0");
        for(j=i; j<dim; ++j)
            mpfr_fprintf(fp," %Rf", self->gmat[self->offset[i] + j], rnd);
        putc('\n', fp);
    }

    fprintf(fp, "MPFR: z:");
    for(i=0; i<dim; ++i)
        mpfr_fprintf(fp, " %Rf", self->z[i]);
    putc('\n', fp);

    fprintf(fp, "MPFR: H:\n");
    for(i=0; i<dim; ++i) {
        for(j=0; j<i; ++j)
            fprintf(fp, " 0");
        for(j=i; j<dim; ++j)
            mpfr_fprintf(fp," %Rf", self->hmat[self->offset[i] + j], rnd);
        putc('\n', fp);
    }
}
