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
static inline int Dbl_equals_allowNonfinite(double x, double y);
int         compareLongs(const void *void_x, const void *void_y);
int         compareDoubles(const void *void_x, const void *void_y);
int         getNumCores(void);
void        dostacktrace(const char *file, int line, FILE * ofp);
double      perturb_ratio_w(double x, double w, gsl_rng * rng);
double      perturb_ratio(double x, gsl_rng * rng);
long double perturb_interval(long double x, long double lo, long double hi,
                             gsl_rng * rng);
int         removeZeroes(int dim, unsigned x[dim]);
void        eprintf(const char *fmt, ...);
FILE       *efopen(const char *restrict name, const char *restrict mode);
void        printBranchTab(double tab[3][3], FILE * fp);
void        unitTstResult(const char *facility, const char *result);
void        tellopt(const char *opt, const char *description);
unsigned    Dbl_first_geq(double val, unsigned len, double v[len]);
long        long_first_geq(long val, long *v, long len);
long        long_last_leq(long val, long *v, long len);
int         strCountSetChunks(const char *str, const char *sep);
void       *memdup(const void *p, size_t n);
double      KLdiverg(int n, const double o[n], const double e[n]);
double      sum_double(int n, const double x[n]);
double      reflect(double x, double lo, double hi);
char       *strlowercase(char *s);
unsigned    strhash(const char *ss);
int         stripchr(char *s, int c);
char       *stripWhiteSpace(char *buff);
char       *stripInternalWhiteSpace(char *buff);
char       *nextWhitesepToken(char **str);
int         tokenize(int dim, char *token[dim], char *s, const char *delim);
void        strReplaceChr(char *s, int a, int b);
double      parseDbl(char *token);
char       *strltrunc(char *s, int n);
void        hdr(const char *msg);
char       *strcenter(const char *text, unsigned width,
                      char *buff, size_t buffsize);
int         readline(int dim, char buff[dim], FILE *fp);

static inline double survival(double t, double twoN);

#  define ERR(code, msg) do{                        \
        char buff[50];                              \
        strerror_r((code), buff, sizeof(buff));     \
        fprintf(stderr,"%s:%s:%d: %s: %d (%s)\n",   \
                __FILE__,__func__,__LINE__,         \
                (msg), (code), buff);               \
        exit(1);                                    \
    }while(0)

#  define DIE(msg) do{                              \
        fprintf(stderr,"%s:%s:%d: %s\n",            \
                __FILE__,__func__,__LINE__, (msg)); \
        exit(EXIT_FAILURE);                         \
    }while(0)

#  define   CHECKMEM(x) do {                                \
        if((x)==NULL) {                                     \
            fprintf(stderr, "%s:%s:%d: NULL pointer\n",     \
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

// If PTR!=NULL, add offset OSET to pointer PTR (if SIGN>0) or
// subtract it (otherwise). Units are sizeof(char) rather than the
// size of the object to which PTR refers.
#define SHIFT_PTR(PTR,OSET,SIGN) do{                                    \
        if((PTR) != NULL) {                                             \
            if((SIGN) > 0)                                              \
                (PTR) = (void *) (((size_t) (PTR)) + ((size_t) (OSET))); \
            else                                                        \
                (PTR) = (void *) (((size_t) (PTR)) - ((size_t) (OSET))); \
        }                                                               \
    }while(0);

/// Return 1 if the relative difference between x and y is less than or
/// equal to 8*DBL_EPSILON.
static inline int Dbl_near(double x, double y) {
    return fabs(x - y) <= fmax(fabs(x), fabs(y)) * 8.0 * DBL_EPSILON;
}

/// Return 1 if x==y, 0 otherwise. Unlike the standard comparison,
/// this one returns 1 if both values are nan, if both are positive infinity,
/// or if both are negative infinity.
static inline int Dbl_equals_allowNonfinite(double x, double y) {
    if(x == y)
        return 1;

    int xtype = fpclassify(x);
    int ytype = fpclassify(y);

    if(xtype==FP_NAN && ytype==FP_NAN)
        return 1;
    if(xtype==FP_INFINITE && ytype==FP_INFINITE) {
        if(x>0.0 && y>0.0)
            return 1;
        if(x<0.0 && y<0.0)
            return 1;
    }
    return 0;
}

/// survival fuction
static inline double survival(double t, double twoN) {
    assert(t >= 0.0);
    assert(twoN > 0.0);
    return exp(-t / twoN);
}

#endif
