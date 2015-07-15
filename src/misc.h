#ifndef GPTREE_MISC_H
#define GPTREE_MISC_H

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_rng.h>

static inline int Dbl_near(double x, double y);
int getNumCores(void);
void dostacktrace(const char *file, int line, FILE * ofp);
double perturb_ratio_w(double x, double w, gsl_rng *rng);
double perturb_ratio(double x, gsl_rng *rng);
long double perturb_interval(long double x, long double lo, long double hi, gsl_rng *rng);
void eprintf(const char *fmt, ...);
void checkmem( /*@null@ */ void *obj, const char *file, int line);
void printBranchTab(double tab[3][3], FILE *fp);
void normItab(double tab[3][3], unsigned nreps);
void printItab(double tab[3][3], FILE *fp);
void printEtab(double x10, double x20, double x21, FILE *fp);
void pictex(double *x, double *y, unsigned n, const char *xlbl,
            const char *ylbl, const char *plotheading, const char *fname);
void pictex_par(double *x, double *y, unsigned n, const char *xlbl,
            const char *ylbl, const char *plotheading, const char *fname);
unsigned Dbl_first_geq(double val, unsigned len, double v[len]);


#define ERR(code, msg) do{\
    fprintf(stderr,"%s:%s:%d: %s %d (%s)\n",\
            __FILE__,__func__,__LINE__,\
            (msg), (code), strerror((code)));   \
    exit(1);\
}while(0)

#define REQUIRE(x,file,lineno) do { \
  if (!(x)) { \
    dostacktrace(__FILE__,__LINE__,stderr); \
    eprintf("ERR@%s:%d->%s:%d: Sanity check FAIL\n",\
            (file),(lineno),__FILE__,__LINE__); \
   }\
} while(0)

/**
 * Return 1 if the relative difference between x and y is less than or
 * equal to DBL_EPSILON.
 */
static inline int Dbl_near(double x, double y) {
    return fabs(x - y) <= fmax(fabs(x), fabs(y)) * DBL_EPSILON;
}

#endif
