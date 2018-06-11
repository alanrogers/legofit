#include "hessian.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double quadsurf(int dim, double dat[dim], int offset[dim], int nterms,
                double tru[nterms], gsl_rng *rng, double sigma);

/**
 * Return the ordinate of a quadratic surface corresponding to the
 * abscissa in array "dat", and the regression parameters in array
 * "tru".
 */
double quadsurf(int dim, double dat[dim], int offset[dim], int nterms,
                double tru[nterms], gsl_rng *rng, double sigma) {
    double b, x, y, z=tru[0];
    int i, j;
    for(i=0; i < dim; ++i)
        z += dat[i]*tru[i+1];
    for(i=0; i<dim; ++i) {
        x = dat[i];
        for(j=i; j<dim; ++j) {
            y = dat[j];
            b = tru[1+dim + offset[j] + i];
            z += b*x*y;
        }
    }
    return z + gsl_ran_gaussian(rng, sigma);
}

int main(void) {
    int i,j,k, dim = 20, nrows=1000;
    int nterms = 1 + dim + (dim*(dim+1))/2;
    const char *fname = "xhessian.tmp";
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, time(NULL));
    double sigma=0.05;

    double tru[nterms];
    int offset[dim];
    tru[0] = 0.0;
    for(i=0; i<dim; ++i) {
        offset[i] = (i*(i+1))/2;
        tru[1+i] = 0.0;
    }

    for(j=0; j<dim; ++j) {
        for(k=j; k<dim; ++k)
            tru[1+dim + offset[k]+j] = -1.0/sqrt(1.0+j+k);
    }

    for(i=0; i<nterms; ++i)
        printf("tru[%2d]=%lg\n", i, tru[i]);

    // Make data file
    double dat[dim];
    FILE *fp = fopen(fname, "w");
    if(fp==NULL) {
        fprintf(stderr,"%s:%d: can't read %s\n",__FILE__,__LINE__,
                fname);
        exit(EXIT_FAILURE);
    }
    fprintf(fp,"%d %d\n", nrows, dim);
    fprintf(fp, "%s", "lnL");
    for(j=0; j<dim; ++j)
        fprintf(fp," theta%d", j);
    putc('\n', fp);
    for(i=0; i<nrows; ++i) {
        if(i==0) {
            for(j=0; j<dim; ++j)
                dat[j] = 0.0;
        }else{
            for(j=0; j<dim; ++j)
                dat[j] = gsl_ran_flat(rng, -1.0, 1.0);
        }
        fprintf(fp, "%0.18g", quadsurf(dim, dat, offset, nterms,
                                       tru, rng, sigma));
        for(j=0; j<dim; ++j)
            fprintf(fp, " %0.18g", dat[j]);
        putc('\n', fp);
    }
    fclose(fp);

    // Construct Hessian matrix
    Hessian hesobj = hessian(fname);
    gsl_matrix *H = hesobj.hessian;
    char **Hparname = hesobj.parname;
    double lnL = hesobj.lnL;

    printf("lnL=%lg\n", lnL);

    double x, ex, err;
    for(i=0; i<dim; ++i) {
        for(j=i; j<dim; ++j) {
            x = gsl_matrix_get(H, i, j);
            ex = tru[1+dim+offset[j]+i];
            if(i==j)
                x /= 2.0;
            err = x - ex;
            printf("[%d,%d]: x=%0.5lf Ex=%0.5lf err=%lg\n",
                   i, j, x, ex, err);
        }
    }

    gsl_matrix_free(H);
    for(i=0; i<dim; ++i)
        free(Hparname[i]);
    free(Hparname);
    unlink(fname);
    return 0;
}
