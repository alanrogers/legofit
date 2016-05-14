#ifndef GPTREE_MISC_H
#  define GPTREE_MISC_H

#  include "typedefs.h"
#  include "binary.h"
#  include "lblndx.h"
#  include <stdio.h>
#  include <math.h>
#  include <float.h>
#  include <assert.h>
#  include <gsl/gsl_rng.h>

static inline int Dbl_near(double x, double y);
int         getNumCores(void);
void        dostacktrace(const char *file, int line, FILE * ofp);
double      perturb_ratio_w(double x, double w, gsl_rng * rng);
double      perturb_ratio(double x, gsl_rng * rng);
long double perturb_interval(long double x, long double lo, long double hi,
                             gsl_rng * rng);
void        eprintf(const char *fmt, ...);
void        checkmem( /*@null@ */ void *obj, const char *file, int line);
void        printBranchTab(double tab[3][3], FILE * fp);
void        normItab(double tab[3][3], unsigned nreps);
void        printItab(double tab[3][3], FILE * fp);
void        printEtab(double x10, double x20, double x21, FILE * fp);
void        pictex(double *x, double *y, unsigned n, const char *xlbl,
                   const char *ylbl, const char *plotheading,
                   const char *fname);
void        pictex_par(double *x, double *y, unsigned n, const char *xlbl,
                       const char *ylbl, const char *plotheading,
                       const char *fname);
void        unitTstResult(const char *facility, const char *result);
void        tellopt(const char *opt, const char *description);
unsigned    Dbl_first_geq(double val, unsigned len, double v[len]);
int         strCountSetChunks(const char *str, const char *sep);
void       *memdup(const void *p, size_t n);
char       *patLbl(size_t n, char buff[n], tipId_t tid, const LblNdx * lblndx);
void        orderpat(int n, unsigned order[n], tipId_t tid[n]);
double      KLdiverg(int n, const double o[n], const double e[n]);
double      sum_double(int n, const double x[n]);

static inline double survival(double t, double twoN);

#  define ERR(code, msg) do{\
    fprintf(stderr,"%s:%s:%d: %s %d (%s)\n",\
            __FILE__,__func__,__LINE__,\
            (msg), (code), strerror((code)));   \
    exit(1);\
}while(0)

#  define   CHECKMEM(x) do {                                  \
        if((x)==NULL) {                                     \
            fprintf(stderr, "%s:%s:%d: allocation error\n", \
                    __FILE__,__func__,__LINE__);            \
            exit(EXIT_FAILURE);                             \
        }                                                   \
    } while(0)

#  define REQUIRE(x,file,lineno) do { \
  if (!(x)) { \
    dostacktrace(__FILE__,__LINE__,stderr); \
    eprintf("ERR@%s:%d->%s:%d: Sanity check FAIL\n",\
            (file),(lineno),__FILE__,__LINE__); \
   }\
} while(0)

/// Return 1 if the relative difference between x and y is less than or
/// equal to DBL_EPSILON.
static inline int Dbl_near(double x, double y) {
    return fabs(x - y) <= fmax(fabs(x), fabs(y)) * DBL_EPSILON;
}

/// survival fuction
static inline double survival(double t, double twoN) {
    assert(t >= 0.0);
    assert(twoN > 0.0);
    return exp(-t / twoN);
}

#endif
