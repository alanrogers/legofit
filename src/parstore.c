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
#include "strparmap.h"
#include "addrparmap.h"
#include "dtnorm.h"
#include "tinyexpr.h"
#include "ptrset.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>

// Set value of constrained variable. Abort with an error message if
// result is NaN.
#define SET_CONSTR(par) do {                                            \
        (par)->value = te_eval((par)->constr);                          \
        if(isnan((par)->value) {                                        \
            fprintf(stderr,"%s:%d: constraint returned NaN\n",          \
                    __FILE__,__LINE__);                                 \
            fprintf(stderr,"formula: %s = %s\n",                        \
                    (par)->name, (par)->formula);                       \
            exit(EXIT_FAILURE);                                         \
        }                                                               \
    }while(0)

#define NEWCODE
#ifdef NEWCODE

struct ParStore {
    int nFree, nUnfree;
    Param freePar[MAXPAR];
    Param unfreePar[MAXPAR];
    StrParMap *byname;               // look up by name
    AddrParMap *byaddr;              // look up by address of value
    te_variable *te_pars;            // linked list
};

/// Set vector of free parameters.
void ParStore_setFreeParams(ParStore *self, int n, double x[n]) {
    assert(n == self->nFree);
    for(int i=0; i < n; ++i)
        self->freePar[i].value = x[i];
}

/// Get vector of free parameters.
void ParStore_getFreeParams(ParStore *self, int n, double x[n]) {
    assert(n == self->nFree);
    for(int i=0; i < n; ++i)
        x[i] = self->freePar[i].value;
}

/// Print a ParStore
void ParStore_print(ParStore *self, FILE *fp) {
    int i;
    fprintf(fp, "Fixed:\n");
    for(i=0; i < self->nUnfree; ++i)
        Param_print(self->unfreePar + i, FIXED, fp);

    ParStore_printFree(self, fp);
    ParStore_printConstrained(self, fp);
}

/// Print free parameter values
void ParStore_printFree(ParStore *self, FILE *fp) {
    int i;
    fprintf(fp, "%5d free:\n", self->nFree);
    for(i=0; i < self->nFree; ++i)
        Param_print(self->freePar+i, FREE, fp);
}

/// Print constrained parameter values
void ParStore_printConstrained(ParStore *self, FILE *fp) {
    int i;
    fprintf(fp, "constrained:\n");
    for(i=0; i < self->nUnfree; ++i)
        Param_print(self->unfreePar+i, CONSTRAINED, fp);
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
    int status;
    ParStore   *new = memdup(old, sizeof(ParStore));
    new->byname = new->byaddr = NULL;

    int i, j=0, k=0;
    for(i = 0; i < new->nFree; ++i) {
        new->byname = StrParMap_insert(new->byname, new->freePar + i);
        new->byaddr = AddrParMap_insert(new->byaddr, new->freePar + i);
        new->te_pars = te_variable_push(new->te_pars, new->freePar[i].name,
                                        &new->freePar[i].value);
    }

    for(i = 0; i < new->nUnfree; ++i) {
        Param *par = new->unfreePar + i;

        // skip parameters that are not CONSTRAINED
        if( par->type != Constrained )
            continue;

        new->te_pars = te_variable_push(new->te_pars, par->name, &par->value);
        par->formula = strdup(old->unfreePar[i].formula);
        par->constr = te_compile(par->formula, new->te_pars, &status);
        if(par->constr == NULL) {
            fprintf(stderr,"%s:%d: parse error\n", __FILE__,__LINE__);
            fprintf(stderr,"  %s\n", par->formula);
            fprintf(stderr,"  %*s^\nError near here\n", status-1, "");
            exit(EXIT_FAILURE);
        }
    }
    ParStore_sanityCheck(new, __FILE__, __LINE__);
    return new;
}

/// Destructor
void ParStore_free(ParStore * self) {
    int         i;

    StrParMap_free(self->byname);
    AddrParMap_free(self->byaddr);

    for(i = 0; i < self->nFree; ++i)
        Param_freePtrs(self->freeParam + i);

    for(i = 0; i < self->nUnfree; ++i)
        Param_freePtrs(self->unfreeParam + i);
    free(self);
}

/// Add a free parameter to ParStore.
void ParStore_addFreePar(ParStore * self, double value,
                         double lo, double hi, const char *name) {
    Param *par = StrParMap_search(self->byname, name);
    if(par) {
        fprintf(stderr,"%s:%d: Duplicate definition of parameter \"%s\".\n",
                __FILE__,__LINE__, name);
        exit(EXIT_FAILURE);
    }

    int i = self->nFree;
    if(++self->nFree >= MAXPAR) {
        fprintf(stderr, "%s:%s:%d: buffer overflow."
                " nFree=%d. MAXPAR=%d."
                " Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__, self->nFree, MAXPAR);
        exit(1);
    }

    Param_init(self->freePar + i, name, value, low, high, FREE);
    self->byname = StrParMap_insert(self->byname, self->freePar + i);
    self->addr = AddrParMap_insert(self->addr, self->freePar + i);

    self->te_pars = te_variable_push(self->te_pars,
                                     self->freePar[i].name,
                                     &self->freePar[i].value);
}

/// Add fixed parameter to ParStore.
void ParStore_addFixedPar(ParStore * self, double value, const char *name) {
    Param *par = StrParMap_search(self->byname, name);
    if(par) {
        fprintf(stderr,"%s:%d: Duplicate definition of parameter \"%s\".\n",
                __FILE__,__LINE__, name);
        exit(EXIT_FAILURE);
    }

    int i = self->nUnfree;
    if(++self->nUnfree >= MAXPAR) {
        fprintf(stderr, "%s:%s:%d: buffer overflow."
                " nUnfree=%d. MAXPAR=%d."
                " Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__, self->nUnfree, MAXPAR);
        exit(EXIT_FAILURE);
    }

    // For fixed parameters, the lower and upper bounds equal the
    // parameter value.
    Param_init(self->unfreePar + i, name, value, value, value, FIXED);
    self->byname = StrParMap_insert(self->byname, self->unfreePar + i);
    self->addr = AddrParMap_insert(self->addr, self->unfreePar + i);
}

/// Add constrained parameter to ParStore.
void ParStore_addConstrainedPar(ParStore * self, const char *str,
                                const char *name) {
    Param *par = StrParMap_search(self->byname, name);
    if(par) {
        fprintf(stderr,"%s:%d: Duplicate definition of parameter \"%s\".\n",
                __FILE__,__LINE__, name);
        exit(EXIT_FAILURE);
    }

    int i = self->nUnfree;
    if(++self->nUnfree >= MAXPAR) {
        fprintf(stderr, "%s:%s:%d: buffer overflow."
                " nUnfree=%d. MAXPAR=%d."
                " Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__, self->nUnfree, MAXPAR);
        exit(EXIT_FAILURE);
    }

    par = self->unfreePar + i;
    Param_init(par, name, value, DBL_MIN, DBL_MAX, CONSTRAINED);
    par->formula = strdup(str);
    CHECKMEM(par->formula);

    self->byname = StrParMap_insert(self->byname, par);
    self->addr = AddrParMap_insert(self->addr, par);

    int status;
    par->constr = te_compile(str, self->te_pars, &status);
    if(self->constr[j] == NULL) {
        fprintf(stderr,"%s:%d: parse error\n", __FILE__,__LINE__);
        fprintf(stderr,"  %s\n", str);
        fprintf(stderr,"  %*s^\nError near here\n", status-1, "");
        exit(EXIT_FAILURE);
    }
    SET_CONSTR(par);
    self->te_pars = te_variable_push(self->te_pars, par->name, &par->value);
}

/// Return the number of free parameters
int ParStore_nFree(ParStore * self) {
    assert(self);
    return self->nFree;
}

/// Return the number of fixed parameters
int ParStore_nFixed(ParStore * self) {
    assert(self);
    int n=0;
    for(int i=0; i < self->nUnfree; ++i)
        if(self->unfreePar[i].type == Fixed)
            ++n;
    return n;
}

/// Return the number of constrained parameters
int ParStore_nConstrained(ParStore * self) {
    assert(self);
    int n=0;
    for(int i=0; i < self->nUnfree; ++i)
        if(self->unfreePar[i].type == Constrained)
            ++n;
    return n;
}

/// Get name of i'th free parameter
const char *ParStore_getNameFree(ParStore * self, int i) {
    assert(i < self->nFree);
    return self->freePar[i].name;
}

/// Return pointer associated with parameter name.
double *ParStore_findPtr(ParStore * self, ParamType *type, const char *name) {
    Param *par = StrParMap_search(self->byname, name);
    if(par == NULL)
        return NULL;
    *type = par->type;
    return &par->value;
}

/// Return 1 if ptr is the address of a constrained parameter; 0
/// otherwise.  
int ParStore_isConstrained(const ParStore *self, const double *ptr) {
    assert(self);
    Param *par = AddrParMap(self->byaddr, ptr);
    if(par && (par->type == Constrained))
        return 1;
    return 0;
}

void ParStore_chkDependencies(ParStore *self, double *valptr, PtrSet *seen) {
    assert(self);

    Param *par = AddrParMap(self->byaddr, valptr);
    
    // Nothing to be done unless par is constrained.
    if( par->type != Constrained )
        return;

    // index of par in array
    int ipar = par - self->constrainedVal;

    // Get list of pointers to parameters on which par depends.
    int len = 100;
    const double * dep[len];
    len = te_dependencies(par->constr, len, dep);

    // Check that each dependendant constrained parameter is in
    // "seen". This implies that dependencies will be set before par
    // is set.
    for(int j=0; j<len; ++j) {
        // dpar is Param structure associated with value pointer dep[j]
        Param *dpar = AddrParMap(self->byaddr, dep[j]);

        if(dpar == NULL) {
            fprintf(stderr,"%s:%d: can't find parameter with address %p\n",
                    __FILE__,__LINE__, dep[j]);
            exit(EXIT_FAILURE);
        }

        if(!(dpar & CONSTRAINED)) {
            // dep[j] not constrained: no problem
            continue;
        }

        // Does par depend on itself?
        if(par == dpar) {
            fprintf(stderr, "%s:%d: Error: \"%s\" depends on itself\n",
                    __FILE__,__LINE__, par->name);
            exit(EXIT_FAILURE);
        }

        // If x and y are constrained parameters and y = f(x), then x
        // must come before y in the .lgo file.
        if(par < dpar) {
            fprintf(stderr, "%s:%d: Error: \"%s\" depends on"
                    " \"%s\" and must come later in\n"
                    " the .lgo file.\n",
                    __FILE__,__LINE__, par->name, dpar->name);
            exit(EXIT_FAILURE);
        }

        if(PtrSet_exists(seen, dep[j])) {
            // dep[i] is set before par: no problem. 
            continue;
        }

        // Pointer dep[j] is not in "seen", which means that it will
        // not be set before "par" is evaluated. Print an error
        // message and abort.
        fprintf(stderr,"%s:%d: Error: \"%s\" depends on \"%s\".\n",
                __FILE__,__LINE__, par->name, dpar->name);
        fprintf(stderr, "   When one constrained parameter depends on"
                " another, the dependent\n"
                "   parameter must appear in the tree, either earlier"
                " on the same\n"
                "   line of descent or on a different line of descent.\n");
        exit(EXIT_FAILURE);
    }
}

void ParStore_sanityCheck(ParStore *self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self, file, line);
    REQUIRE(self->nFree >= 0, file, line);
    REQUIRE(self->nUnfree >= 0, file, line);
    REQUIRE(self->nFree < MAXPAR, file, line);
    REQUIRE(self->nUnfree < MAXPAR, file, line);

    // For each name: (1) make sure it's a legal name;
    // (2) get the pointer associated with that name,
    int i;
    char *s;
    double *ptr;
    ParamType ptype;
    Param *par;
    for(i=0; i < self->nFree; ++i) {
        par = self->freePar + i;
        REQUIRE(NULL != par->name, file, line);
        REQUIRE(legalName(par->name), file, line);
        ptr = ParKeyVal_get(self->pkv, &ptype, s);
        REQUIRE(ptr != NULL, file, line);
        REQUIRE(ptype == Fixed, file, line);
        REQUIRE(ptr == self->fixedVal+i, file, line);
    }
    for(i=0; i < self->nFree; ++i) {
        s = self->nameFree[i];
        REQUIRE(NULL != s, file, line);
        REQUIRE(legalName(s), file, line);
        ptr = ParKeyVal_get(self->pkv, &ptype, s);
        REQUIRE(ptr != NULL, file, line);
        REQUIRE(ptype == Free, file, line);
        REQUIRE(ptr == self->freeVal+i, file, line);
    }
    for(i=0; i < self->nGaussian; ++i) {
        s = self->nameGaussian[i];
        REQUIRE(NULL != s, file, line);
        REQUIRE(legalName(s), file, line);
        ptr = ParKeyVal_get(self->pkv, &ptype, s);
        REQUIRE(ptr != NULL, file, line);
        REQUIRE(ptype == Gaussian, file, line);
        REQUIRE(ptr == self->gaussianVal+i, file, line);
    }
    for(i=0; i < self->nConstrained; ++i) {
        s = self->nameConstrained[i];
        REQUIRE(NULL != s, file, line);
        REQUIRE(legalName(s), file, line);
        ptr = ParKeyVal_get(self->pkv, &ptype, s);
        REQUIRE(ptr != NULL, file, line);
        REQUIRE(ptype == Constrained, file, line);
        REQUIRE(ptr == self->constrainedVal+i, file, line);
    }
    ParKeyVal_sanityCheck(self->pkv, file, line);
#endif
}

#else

/// Return 1 if two ParStore objects are equal; 0 otherwise.
int         ParStore_equals(ParStore *lhs, ParStore *rhs) {
    if(lhs == rhs)
        return 1;
    if(lhs->nFixed != rhs->nFixed) {
        return 0;
    }
    if(lhs->nFree != rhs->nFree) {
        return 0;
    }
    if(lhs->nGaussian != rhs->nGaussian) {
        return 0;
    }
    if(lhs->nConstrained != rhs->nConstrained) {
        return 0;
    }
    if(0 != memcmp(lhs->fixedVal, rhs->fixedVal,
                   lhs->nFixed*sizeof(lhs->fixedVal[0]))) {
        return 0;
    }
    if(0 != memcmp(lhs->freeVal, rhs->freeVal,
                   lhs->nFree*sizeof(lhs->freeVal[0]))) {
        return 0;
    }
    if(ParStore_constrain(lhs)) {
        fprintf(stderr,"%s:%d: free parameters violate constraints\n",
                __FILE__,__LINE__);
    }
    if(ParStore_constrain(rhs)) {
        fprintf(stderr,"%s:%d: free parameters violate constraints\n",
                __FILE__,__LINE__);
    }
    if(0 != memcmp(lhs->constrainedVal, rhs->constrainedVal,
                   lhs->nConstrained*sizeof(lhs->constrainedVal[0]))) {
        return 0;
    }
    if(0 != memcmp(lhs->loFree, rhs->loFree,
                   lhs->nFree*sizeof(lhs->loFree[0]))) {
        return 0;
    }
    if(0 != memcmp(lhs->hiFree, rhs->hiFree,
                   lhs->nFree*sizeof(lhs->hiFree[0]))) {
        return 0;
    }
    if(0 != memcmp(lhs->mean, rhs->mean,
                   lhs->nGaussian*sizeof(lhs->mean[0]))) {
        return 0;
    }
    if(0 != memcmp(lhs->sd, rhs->sd,
                   lhs->nGaussian*sizeof(lhs->sd[0]))) {
        return 0;
    }
    int i;
    for(i=0; i < lhs->nFixed; ++i)
        if(0 != strcmp(lhs->nameFixed[i], rhs->nameFixed[i])){
            return 0;
        }
    for(i=0; i < lhs->nFree; ++i)
        if(0 != strcmp(lhs->nameFree[i], rhs->nameFree[i])) {
            return 0;
        }
    for(i=0; i < lhs->nGaussian; ++i)
        if(0 != strcmp(lhs->nameGaussian[i], rhs->nameGaussian[i])) {
            return 0;
        }
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
    SET_CONSTR(i);
    assert(*ptr == self->constrainedVal[i]);
}

/// First check to see that free parameters obey boundary constraints.
/// If not, then return 1. Otherwise, set values of all constrained
/// parameters and return 0.
int ParStore_constrain(ParStore *self) {
    int i;
    for(i=0; i < self->nFree; ++i) {
        if(self->freeVal[i] < self->loFree[i]
           || self->freeVal[i] > self->hiFree[i])
            return 1;
    }
    for(i=0; i < self->nConstrained; ++i) {
        SET_CONSTR(i);
    }
    return 0;
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

#endif
