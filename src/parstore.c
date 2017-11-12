/**
 * @file parstore.c
 * @brief Manage a vector of parameters
 *
 * This class solves the following problem. One algorithm requires
 * that parameters be distributed throughout a network of
 * nodes. Another requires that they all be collected into a
 * vector. To accomplish both goals, ParStore maintains a vector of
 * parameter values, together with high and low bounds on those
 * values. When a new parameter is added to the ParStore, a pointer to
 * that value is returned, so that it can be stored in the distributed
 * data structure.
 *
 * In fact, it maintains 3 such vectors: one each for free parameters,
 * fixed parameters, and Gaussian parameters. ParStore knows how to
 * perturb Gaussian parameters by sampling from a truncated Gaussian
 * distribution. Fixed parameters never change. Free ones are manipulated
 * from outside via function calls.
 *
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "parstore.h"
#include "parkeyval.h"
#include "dtnorm.h"
#include "tinyexpr.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

struct ParStore {
    int         nFixed, nFree, nGaussian, nConstrained; // num pars
    double      loFree[MAXPAR]; // lower bounds
    double      hiFree[MAXPAR]; // upper bounds
    char       *nameFixed[MAXPAR];  // Parameter names
    char       *nameFree[MAXPAR];   // Parameter names
    char       *nameGaussian[MAXPAR];    // Parameter names
    char       *nameConstrained[MAXPAR]; // Parameter names
    ParKeyVal  *pkv;           // linked list of name/ptr pairs
    double      fixedVal[MAXPAR];   // parameter values
    double      freeVal[MAXPAR];    // parameter values
    double      gaussianVal[MAXPAR]; // parameter values
    double      constrainedVal[MAXPAR]; // parameter values
    double      mean[MAXPAR];        // Gaussian means
    double      sd[MAXPAR];          // Gaussian standard deviations
    te_expr    *constr[MAXPAR];      // controls constrainedVal entries
    te_variables te_pars[MAXPAR];
    char       *formulas[MAXPAR];    // formulas of constrained vars
};

static int compareDblPtrs(const void *void_x, const void *void_y);
static int compareDbls(const void *void_x, const void *void_y);
static inline int chrcount(const char *s, char c);

/// Count the number of copies of character c in string s
static inline int chrcount(const char *s, char c) {
    if(s==NULL)
        return 0;
    int n=0;
    while(*s != '\0') {
        if(*s++ == c)
            ++n;
    }
    return n;
}

/// Return <0, 0, or >0, as x is <, ==, or > y.
static int compareDblPtrs(const void *void_x, const void *void_y) {
    double * const * x = (double * const *) void_x;
    double * const * y = (double * const *) void_y;
    return **x - **y;
}

/// Return <0, 0, or >0, as x is <, ==, or > y.
static int compareDbls(const void *void_x, const void *void_y) {
    const double * x = (const double *) void_x;
    const double * y = (const double *) void_y;
    return *x - *y;
}

/// Set vector of free parameters.
void ParStore_setFreeParams(ParStore *self, int n, double x[n]) {
    assert(n == self->nFree);
    memcpy(self->freeVal, x, n*sizeof(double));
}

/// Get vector of free parameters.
void ParStore_getFreeParams(ParStore *self, int n, double x[n]) {
    assert(n == self->nFree);
    memcpy(x, self->freeVal, n*sizeof(double));
}

/// Print a ParStore
void ParStore_print(ParStore *self, FILE *fp) {
    int i;
    fprintf(fp, "%5d fixed:\n", self->nFixed);
    for(i=0; i < self->nFixed; ++i)
        fprintf(fp, "   %8s = %lg\n", self->nameFixed[i], self->fixedVal[i]);
    fprintf(fp, "%5d Gaussian:\n", self->nGaussian);
    for(i=0; i < self->nGaussian; ++i)
        fprintf(fp, "   %8s = Gaussian(%lg, %lg)\n", self->nameGaussian[i],
                self->mean[i], self->sd[i]);
    ParStore_printFree(self, fp);
    ParStore_printConstrained(self, fp);
}

/// Print free parameter values
void ParStore_printFree(ParStore *self, FILE *fp) {
    int i;
    fprintf(fp, "%5d free:\n", self->nFree);
    for(i=0; i < self->nFree; ++i)
        fprintf(fp, "   %8s = %lg\n", self->nameFree[i], self->freeVal[i]);
}

/// Print constrained parameter values
void ParStore_printConstrained(ParStore *self, FILE *fp) {
    int i;
    fprintf(fp, "%5d constrained:\n", self->nConstrained);
    for(i=0; i < self->nConstrained; ++i) {
        fprintf(fp, "   %8s = %lg = ", self->nameConstrained[i],
                self->constrainedVal[i]);
        te_print(self->constr[i], fp);
    }
}

/// Constructor
ParStore *ParStore_new(void) {
    ParStore   *self = malloc(sizeof(ParStore));
    CHECKMEM(self);
    memset(self, 0, sizeof(ParStore));
    ParStore_sanityCheck(self, __FILE__, __LINE__);
    return self;
}

/// Duplicate a ParStore
ParStore *ParStore_dup(const ParStore * old) {
    assert(old);
    ParStore   *new = memdup(old, sizeof(ParStore));
    new->pkv = NULL;

    int         i;
    for(i = 0; i < new->nFree; ++i) {
        new->nameFree[i] = strdup(old->nameFree[i]);
        new->pkv = ParKeyVal_add(new->pkv, new->nameFree[i],
                                  new->freeVal + i, Free);
        new->te_pars[i].name = new->nameFree[i];
        new->te_pars[i].address = new->freeVal + i;
    }

    for(i = 0; i < new->nFixed; ++i) {
        new->nameFixed[i] = strdup(old->nameFixed[i]);
        new->pkv = ParKeyVal_add(new->pkv, new->nameFixed[i],
                                  new->fixedVal + i, Fixed);
    }

    for(i = 0; i < new->nGaussian; ++i) {
        new->nameGaussian[i] = strdup(old->nameGaussian[i]);
        new->pkv = ParKeyVal_add(new->pkv, new->nameGaussian[i],
                                  new->gaussianVal + i, Gaussian);
    }

    for(i = 0; i < new->nConstrained; ++i) {
        new->nameConstrained[i] = strdup(old->nameConstrained[i]);
        new->pkv = ParKeyVal_add(new->pkv, new->nameConstrained[i],
                                  new->constrainedVal + i, Constrained);
        new->formulas[i] = strdup(old->formulas[i]);
        new->constr[i] = te_compile(new->formulas[i], new->te_pars,
                                    new->nFree, &status);
        if(new->constr[i] == NULL) {
            fprintf(stderr,"%s:%d: parse error\n", __FILE__,__LINE__);
            fprintf(stderr,"  %s\n", new->formulas[i]);
            fprintf(stderr,"  %*s^\nError near here", status-1, "");
            exit(EXIT_FAILURE);
        }
    }
    ParStore_sanityCheck(new, __FILE__, __LINE__);
    return new;
}

/// Destructor
void ParStore_free(ParStore * self) {
    int         i;

    for(i = 0; i < self->nFixed; ++i)
        free(self->nameFixed[i]);

    for(i = 0; i < self->nFree; ++i)
        free(self->nameFree[i]);

    for(i = 0; i < self->nGaussian; ++i)
        free(self->nameGaussian[i]);

    for(i=0; i < self->nConstrained; ++i) {
        free(self->nameConstrained[i]);
        free(self->formulas[i]);
        te_free(self->constr[i]);
    }

    ParKeyVal_free(self->pkv);
    free(self);
}

/// Add free parameter to ParStore.
void ParStore_addFreePar(ParStore * self, double value,
                         double lo, double hi, const char *name) {
    int         i = self->nFree;
    ParamStatus pstat;
    if(NULL != ParKeyVal_get(self->pkv, &pstat, name)) {
        fprintf(stderr,"%s:%d: Duplicate definition of parameter \"%s\".\n",
                __FILE__,__LINE__, name);
        exit(EXIT_FAILURE);
    }

    if(++self->nFree >= MAXPAR) {
        fprintf(stderr, "%s:%s:%d: buffer overflow."
                " nFree=%d. MAXPAR=%d."
                " Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__, self->nFree, MAXPAR);
        exit(1);
    }

    if(value < lo || value > hi) {
        fprintf(stderr,"%s:%s:%d: value (%lf) not in range [%lf,%lf]\n",
                __FILE__, __func__, __LINE__, value, lo, hi);
        exit(1);
    }

    self->freeVal[i] = value;
    self->loFree[i] = lo;
    self->hiFree[i] = hi;
    self->nameFree[i] = strdup(name);
    CHECKMEM(self->nameFree[i]);
    self->te_pars[i].name = self->nameFree[i];
    self->te_pars[i].address = self->freeVal + i;

    // Linked list associates pointer with parameter name.
    self->pkv = ParKeyVal_add(self->pkv, name, self->freeVal + i,
							   Free);
}

/// Add Gaussian parameter to ParStore.
void ParStore_addGaussianPar(ParStore * self, double mean, double sd,
                             const char *name) {
    int         i = self->nGaussian;
    ParamStatus pstat;
    if(NULL != ParKeyVal_get(self->pkv, &pstat, name)) {
        fprintf(stderr,"%s:%d: Duplicate definition of parameter \"%s\".\n",
                __FILE__,__LINE__, name);
        exit(EXIT_FAILURE);
    }

    if(++self->nGaussian >= MAXPAR) {
        fprintf(stderr, "%s:%s:%d: buffer overflow."
                " nGaussian=%d. MAXPAR=%d."
                " Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__, self->nGaussian, MAXPAR);
        exit(1);
    }

    assert(mean >= 0.0);
    assert(sd >= 0.0);

    self->gaussianVal[i] = self->mean[i] = mean;
    self->sd[i] = sd;
    self->nameGaussian[i] = strdup(name);
    CHECKMEM(self->nameGaussian[i]);

    // Linked list associates pointer with parameter name.
    self->pkv = ParKeyVal_add(self->pkv, name, self->gaussianVal + i,
							   Gaussian);
}

/// Add fixed parameter to ParStore.
void ParStore_addFixedPar(ParStore * self, double value, const char *name) {
    int         i = self->nFixed;
    ParamStatus pstat;
    if(NULL != ParKeyVal_get(self->pkv, &pstat, name)) {
        fprintf(stderr,"%s:%d: Duplicate definition of parameter \"%s\".\n",
                __FILE__,__LINE__, name);
        exit(EXIT_FAILURE);
    }

    if(++self->nFixed >= MAXPAR)
        eprintf("%s:%s:%d: buffer overflow."
                " nFixed=%d. MAXPAR=%d."
                " Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__, self->nFixed, MAXPAR);

    self->fixedVal[i] = value;
    self->nameFixed[i] = strdup(name);
    CHECKMEM(self->nameFixed[i]);

    // Linked list associates pointer with parameter name.
    self->pkv = ParKeyVal_add(self->pkv, name, self->fixedVal + i,
		Fixed);
}

/// Add constrained parameter to ParStore.
void ParStore_addConstrainedPar(ParStore * self, char *str,
                                const char *name) {
    int status, i = self->nConstrained;
    ParamStatus pstat;
    if(NULL != ParKeyVal_get(self->pkv, &pstat, name)) {
        fprintf(stderr,"%s:%d: Duplicate definition of parameter \"%s\".\n",
                __FILE__,__LINE__, name);
        exit(EXIT_FAILURE);
    }

    if(++self->nConstrained >= MAXPAR) {
        fprintf(stderr, "%s:%s:%d: buffer overflow."
                " nConstrained=%d. MAXPAR=%d."
                " Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__, self->nConstrained, MAXPAR);
        exit(EXIT_FAILURE);
    }

    self->formulas[i] = strdup(str);
    self->nameConstrained[i] = strdup(name);
    CHECKMEM(self->nameConstrained[i]);

    self->pkv = ParKeyVal_add(self->pkv, name, self->constrainedVal + i,
		Constrained);
    self->constr[i] = te_compile(str, self->te_pars, self->nFree, &status);
    if(self->constr[i] == NULL) {
        fprintf(stderr,"%s:%d: parse error\n", __FILE__,__LINE__);
        fprintf(stderr,"  %s\n", str);
        fprintf(stderr,"  %*s^\nError near here", status-1, "");
        exit(EXIT_FAILURE);
    }
    self->constrainedVal[i] = te_eval(self->constr[i]);
}

/// Return the number of fixed parameters
int ParStore_nFixed(ParStore * self) {
    assert(self);
    return self->nFixed;
}

/// Return the number of free parameters
int ParStore_nFree(ParStore * self) {
    assert(self);
    return self->nFree;
}

/// Return the number of Gaussian parameters
int ParStore_nGaussian(ParStore * self) {
    assert(self);
    return self->nGaussian;
}

/// Return the number of constrained parameters
int ParStore_nConstrained(ParStore * self) {
    assert(self);
    return self->nConstrained;
}

/// Get value of i'th fixed parameter
double ParStore_getFixed(ParStore * self, int i) {
    assert(i < self->nFixed);
    return self->fixedVal[i];
}

/// Get value of i'th free parameter
double ParStore_getFree(ParStore * self, int i) {
    assert(i < self->nFree);
    return self->freeVal[i];
}

/// Get value of i'th Gaussian parameter
double ParStore_getGaussian(ParStore * self, int i) {
    assert(i < self->nGaussian);
    return self->gaussianVal[i];
}

/// Get name of i'th fixed parameter
const char *ParStore_getNameFixed(ParStore * self, int i) {
    assert(i < self->nFixed);
    return self->nameFixed[i];
}

/// Get name of i'th free parameter
const char *ParStore_getNameFree(ParStore * self, int i) {
    assert(i < self->nFree);
    return self->nameFree[i];
}

/// Get name of i'th Gaussian parameter
const char *ParStore_getNameGaussian(ParStore * self, int i) {
    assert(i < self->nGaussian);
    return self->nameGaussian[i];
}

/// Set value of i'th free parameter
void ParStore_setFree(ParStore * self, int i, double value) {
    assert(i < self->nFree);
    self->freeVal[i] = value;
}

/// Return low bound of i'th free parameter
double ParStore_loFree(ParStore * self, int i) {
    assert(i < self->nFree);
    return self->loFree[i];
}

/// Return high bound of i'th free parameter
double ParStore_hiFree(ParStore * self, int i) {
    assert(i < self->nFree);
    return self->hiFree[i];
}

/// Return pointer to array of lower bounds of free parameters
double     *ParStore_loBounds(ParStore * self) {
    return &self->loFree[0];
}

/// Return pointer to array of upper bounds of free parameters
double     *ParStore_upBounds(ParStore * self) {
    return &self->hiFree[0];
}

/// Return pointer associated with parameter name.
double     *ParStore_findPtr(ParStore * self, ParamStatus *pstat,
							 const char *name) {
    return ParKeyVal_get(self->pkv, pstat, name);
}

void ParStore_sanityCheck(ParStore *self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self, file, line);
    REQUIRE(self->nFixed >= 0, file, line);
    REQUIRE(self->nFree >= 0, file, line);
    REQUIRE(self->nGaussian >= 0, file, line);
    REQUIRE(self->nFixed < MAXPAR, file, line);
    REQUIRE(self->nFree < MAXPAR, file, line);
    REQUIRE(self->nGaussian < MAXPAR, file, line);

    // For each name: (1) make sure it's a legal name;
    // (2) get the pointer associated with that name,
    int i;
    char *s;
    double *ptr;
    ParamStatus pstat;
    for(i=0; i < self->nFixed; ++i) {
        s = self->nameFixed[i];
        REQUIRE(NULL != s, file, line);
        REQUIRE(legalName(s), file, line);
        ptr = ParKeyVal_get(self->pkv, &pstat, s);
        REQUIRE(ptr != NULL, file, line);
        REQUIRE(pstat == Fixed, file, line);
        REQUIRE(ptr == self->fixedVal+i, file, line);
    }
    for(i=0; i < self->nFree; ++i) {
        s = self->nameFree[i];
        REQUIRE(NULL != s, file, line);
        REQUIRE(legalName(s), file, line);
        ptr = ParKeyVal_get(self->pkv, &pstat, s);
        REQUIRE(ptr != NULL, file, line);
        REQUIRE(pstat == Free, file, line);
        REQUIRE(ptr == self->freeVal+i, file, line);
    }
    for(i=0; i < self->nGaussian; ++i) {
        s = self->nameGaussian[i];
        REQUIRE(NULL != s, file, line);
        REQUIRE(legalName(s), file, line);
        ptr = ParKeyVal_get(self->pkv, &pstat, s);
        REQUIRE(ptr != NULL, file, line);
        REQUIRE(pstat == Gaussian, file, line);
        REQUIRE(ptr == self->gaussianVal+i, file, line);
    }
    for(i=0; i < self->nConstrained; ++i) {
        s = self->nameConstrained[i];
        REQUIRE(NULL != s, file, line);
        REQUIRE(legalName(s), file, line);
        ptr = ParKeyVal_get(self->pkv, &pstat, s);
        REQUIRE(ptr != NULL, file, line);
        REQUIRE(pstat == Constrained, file, line);
        REQUIRE(ptr == self->constrainedVal+i, file, line);
    }
    ParKeyVal_sanityCheck(self->pkv, file, line);
#endif
}

/// Return 1 if two ParStore objects are equal; 0 otherwise.
int         ParStore_equals(const ParStore *lhs, const ParStore *rhs) {
    if(lhs == rhs)
        return 1;
    if(lhs->nFixed != rhs->nFixed)
        return 0;
    if(lhs->nFree != rhs->nFree)
        return 0;
    if(lhs->nGaussian != rhs->nGaussian)
        return 0;
    if(lhs->nConstrained != rhs->nConstrained)
        return 0;
    if(0 != memcmp(lhs->fixedVal, rhs->fixedVal,
                   lhs->nFixed*sizeof(lhs->fixedVal[0])))
        return 0;
    if(0 != memcmp(lhs->freeVal, rhs->freeVal,
                   lhs->nFree*sizeof(lhs->freeVal[0])))
        return 0;
    if(0 != memcmp(lhs->constrainedVal, rhs->constrainedVal,
                   lhs->nConstrained*sizeof(lhs->constrainedVal[0])))
        return 0;
    if(0 != memcmp(lhs->loFree, rhs->loFree,
                   lhs->nFree*sizeof(lhs->loFree[0])))
        return 0;
    if(0 != memcmp(lhs->hiFree, rhs->hiFree,
                   lhs->nFree*sizeof(lhs->hiFree[0])))
        return 0;
    if(0 != memcmp(lhs->mean, rhs->mean,
                   lhs->nGaussian*sizeof(lhs->mean[0])))
        return 0;
    if(0 != memcmp(lhs->sd, rhs->sd,
                   lhs->nGaussian*sizeof(lhs->sd[0])))
        return 0;
    int i;
    for(i=0; i < lhs->nFixed; ++i)
        if(0 != strcmp(lhs->nameFixed[i], rhs->nameFixed[i]))
            return 0;
    for(i=0; i < lhs->nFree; ++i)
        if(0 != strcmp(lhs->nameFree[i], rhs->nameFree[i]))
            return 0;
    for(i=0; i < lhs->nGaussian; ++i)
        if(0 != strcmp(lhs->nameGaussian[i], rhs->nameGaussian[i]))
            return 0;
    return ParKeyVal_equals(lhs->pkv, rhs->pkv);
}

/// Set Gaussian parameter by sampling from a truncated Gaussian
/// distribution. If ptr points to a Gaussian parameter, a new
/// value of that parameter is drawn from a truncated Gaussian. If the
/// pointer doesn't point to a Gaussian parameter, then the function
/// returns without doing anything.
/// @param[inout] self ParStore object to be modified
/// @param[inout] ptr pointer to parameter to be modified
/// @param[in] low low end of truncation interval
/// @param[in] high high end of truncation interval
/// @param[inout] rng GSL random number generator
void ParStore_sample(ParStore *self, double *ptr, double low, double high,
                    gsl_rng * rng) {
    // If ptr isn't a Gaussian parameter, then return immediately
    if(ptr < self->gaussianVal
       || ptr >= self->gaussianVal + MAXPAR)
        return;

    // index of Gaussian parameter
    unsigned i = ptr - self->gaussianVal;
    assert(i < self->nGaussian);

    // sample from doubly-truncated normal distribution
    self->gaussianVal[i] = dtnorm(self->mean[i], self->sd[i], low, high, rng);
    assert(*ptr == self->gaussianVal[i]);
}

/// If ptr points to a constrained parameter, then set its value.
void ParStore_constrain_ptr(ParStore *self, double *ptr) {
    // If ptr isn't a constrained parameter, then return immediately
    if(ptr < self->constrainedVal
       || ptr >= self->constrainedVal + MAXPAR)
        return;

    // index of constrained parameter
    unsigned i = ptr - self->constrainedVal;
    assert(i < self->nConstrained);

    // set value of constrained parameter
    self->constrainedVal[i] = te_eval(self->constr[i]);
    assert(*ptr == self->constrainedVal[i]);
}

/// Set values of all constrained parameters
void ParStore_constrain(ParStore *self) {
    int i;
    for(i=0; i < self->nConstrained; ++i)
        self->constrainedVal[i] = te_eval(self->constr[i]);
}

/// Make sure Bounds object is sane.
void Bounds_sanityCheck(Bounds * self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self, file, line);
    REQUIRE(self->lo_twoN >= 0.0, file, line);
    REQUIRE(self->lo_twoN < self->hi_twoN, file, line);
    REQUIRE(self->lo_t >= 0.0, file, line);
    REQUIRE(self->lo_t < self->hi_t, file, line);
#endif
}

/// Return 1 if two Bounds objects are equal; 0 otherwise.
int         Bounds_equals(const Bounds *lhs, const Bounds *rhs) {
    if(lhs == rhs)
        return 1;
    return lhs->lo_twoN == rhs->lo_twoN
        && lhs->hi_twoN == rhs->hi_twoN
        && lhs->lo_t == rhs->lo_t
        && lhs->hi_t == rhs->hi_t;
}

#ifdef TEST

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {

    int verbose=0;

    if(argc > 1) {
        if(argc!=2 || 0!=strcmp(argv[1], "-v")) {
            fprintf(stderr,"usage: xterm [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    double x=1.0, y=1.0, z=2.0;
    double *px=&x, *py=&y, *pz=&z;
    double **ppx=&px, **ppy=&py, **ppz=&pz;

    assert(0 == compareDbls(px, py));
    assert(0 == compareDbls(px, px));
    assert(0 > compareDbls(px, pz));
    assert(0 < compareDbls(pz, py));

    unitTstResult("compareDbls", "OK");

    assert(0 == compareDblPtrs(ppx, ppy));
    assert(0 == compareDblPtrs(ppx, ppx));
    assert(0 > compareDblPtrs(ppx, ppz));
    assert(0 < compareDblPtrs(ppz, ppy));

    unitTstResult("compareDblPtrs", "OK");

    ParKeyVal *pkv = NULL;
    pkv = ParKeyVal_add(pkv, "x", &x, Free);
    pkv = ParKeyVal_add(pkv, "y", &y, Free);
    pkv = ParKeyVal_add(pkv, "z", &z, Free);

}
#endif
