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
#include "param.h"
#include "tinyexpr.h"
#include "ptrset.h"
#include "misc.h"
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
        if(isnan((par)->value)) {                                       \
            fprintf(stderr,"%s:%d: constraint returned NaN\n",          \
                    __FILE__,__LINE__);                                 \
            fprintf(stderr,"formula: %s = %s\n",                        \
                    (par)->name, (par)->formula);                       \
            exit(EXIT_FAILURE);                                         \
        }                                                               \
    }while(0)

struct ParStore {
    int nPar;

    // These 3 pointers are the heads of linked lists of free
    // parameters, fixed parameters, and constrained parameters. 
    // The objects in these lists can also be accessed through
    // array "vec".
    Param *freePar, *fixedPar, *constrainedPar;

    // Map names and addresses of values to parameter objects
    StrParMap *byname;          // look up by name
    AddrParMap *byaddr;         // look up by address of value

    te_variable *te_pars;       // for tinyexpr.c

    // This is where the parameter objects are allocated.
    Param vec[MAXPAR];
};

/// Return the number of free parameters
int ParStore_nFree(ParStore * self) {
    assert(self);
    int n = 0;
    for(Param * par = self->freePar; par != NULL; par = par->next)
        ++n;
    return n;
}

/// Return the number of fixed parameters
int ParStore_nFixed(ParStore * self) {
    assert(self);
    int n = 0;
    for(Param * par = self->fixedPar; par != NULL; par = par->next)
        ++n;
    return n;
}

/// Return the number of constrained parameters
int ParStore_nConstrained(ParStore * self) {
    assert(self);
    int n = 0;
    for(Param * par = self->constrainedPar; par != NULL; par = par->next)
        ++n;
    return n;
}

/// Set vector of free parameters.
void ParStore_setFreeParams(ParStore * self, int n, double x[n]) {
    assert(self);
    int i = 0;
    for(Param * par = self->freePar; par != NULL; par = par->next) {
        assert(i < n);
        par->value = x[i++];
    }
    assert(i == n);
}

/// Get vector of free parameters.
void ParStore_getFreeParams(ParStore * self, int n, double x[n]) {
    assert(self);
    int i = 0;
    for(Param * par = self->freePar; par != NULL; par = par->next) {
        assert(i < n);
        x[i++] = par->value;
    }
    assert(i == n);
}

/// Print a ParStore
void ParStore_print(ParStore * self, FILE * fp) {
    fprintf(fp, "Fixed:\n");
    for(Param *par = self->fixedPar; par != NULL; par = par->next)
        Param_print(par, fp);

    ParStore_printFree(self, fp);
    ParStore_printConstrained(self, fp);
}

/// Print free parameter values
void ParStore_printFree(ParStore * self, FILE * fp) {
    fprintf(fp, "Free:\n");
    for(Param *par = self->freePar; par != NULL; par = par->next)
        Param_print(par, fp);
}

/// Print constrained parameter values
void ParStore_printConstrained(ParStore * self, FILE * fp) {
    fprintf(fp, "Constrained:\n");
    for(Param *par = self->constrainedPar; par != NULL; par = par->next)
        Param_print(par, fp);
}

/// Constructor
ParStore *ParStore_new(void) {
    ParStore *self = malloc(sizeof(ParStore));
    CHECKMEM(self);
    memset(self, 0, sizeof(ParStore));
    ParStore_sanityCheck(self, __FILE__, __LINE__);
    return self;
}

/// Duplicate a ParStore
ParStore *ParStore_dup(const ParStore * old) {
    assert(old);
    int status;
    ParStore *new = memdup(old, sizeof(ParStore));
    new->freePar = new->fixedPar = new->constrainedPar = NULL;
    new->byname = NULL;
    new->byaddr = NULL;
    new->te_pars = NULL;

    for(int i = 0; i < new->nPar; ++i) {
        Param *par = new->vec + i;
        const Param *opar = old->vec + i;
        Param_copy(par, opar);
        new->byname = StrParMap_insert(new->byname, par);
        new->byaddr = AddrParMap_insert(new->byaddr, par);
        if(par->type & FREE) {
            new->freePar = Param_push(new->freePar, par);
        }else if(par->type & FIXED) {
            new->fixedPar = Param_push(new->fixedPar, par);
        }else if(par->type & CONSTRAINED) {
            new->constrainedPar = Param_push(new->constrainedPar, par);
            par->constr = te_compile(par->formula, new->te_pars, &status);
            if(par->constr == NULL) {
                fprintf(stderr, "%s:%d: parse error\n", __FILE__, __LINE__);
                fprintf(stderr, "  %s\n", par->formula);
                fprintf(stderr, "  %*s^\nError near here\n", status - 1, "");
                exit(EXIT_FAILURE);
            }
        }else{
            DIE("Illegal type field within ParStore");
        }
        new->te_pars =
            te_variable_push(new->te_pars, par->name, &par->value);
    }

    ParStore_sanityCheck(new, __FILE__, __LINE__);
    return new;
}

/// Destructor
void ParStore_free(ParStore * self) {
    assert(self);
    StrParMap_free(self->byname);
    AddrParMap_free(self->byaddr);
    te_variable_free(self->te_pars);

    for(int i = 0; i < self->nPar; ++i)
        Param_freePtrs(self->vec + i);

    free(self);
}

/// Add a free parameter to ParStore.
void ParStore_addFreePar(ParStore * self, double value, double lo,
                         double hi, const char *name, unsigned type) {
    assert(self);
    Param *par = StrParMap_search(self->byname, name);
    if(par) {
        fprintf(stderr, "%s:%d: Duplicate definition of parameter \"%s\".\n",
                __FILE__, __LINE__, name);
        exit(EXIT_FAILURE);
    }

    int i = self->nPar++;
    if(i == MAXPAR) {
        fprintf(stderr, "%s:%s:%d: buffer overflow."
                " nPar=%d. MAXPAR=%d."
                " Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__, self->nPar, MAXPAR);
        exit(1);
    }

    par = self->vec + i;
    Param_init(par, name, value, lo, hi, type|FREE);
    self->byname = StrParMap_insert(self->byname, par);
    self->byaddr = AddrParMap_insert(self->byaddr, par);
    self->freePar = Param_push(self->freePar, par);
    self->te_pars = te_variable_push(self->te_pars, par->name, &par->value);
    ParStore_sanityCheck(self, __FILE__, __LINE__);
}

/// Add fixed parameter to ParStore.
void ParStore_addFixedPar(ParStore * self, double value, const char *name,
                          unsigned type) {
    Param *par = StrParMap_search(self->byname, name);
    if(par) {
        fprintf(stderr, "%s:%d: Duplicate definition of parameter \"%s\".\n",
                __FILE__, __LINE__, name);
        exit(EXIT_FAILURE);
    }

    int i = self->nPar++;
    if(i == MAXPAR) {
        fprintf(stderr, "%s:%s:%d: buffer overflow."
                " nPar=%d. MAXPAR=%d."
                " Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__, self->nPar, MAXPAR);
        exit(1);
    }

    par = self->vec + i;

    // For fixed parameters, the lower and upper bounds equal the
    // parameter value.
    Param_init(par, name, value, value, value, type | FIXED);
    self->byname = StrParMap_insert(self->byname, par);
    self->byaddr = AddrParMap_insert(self->byaddr, par);
    self->fixedPar = Param_push(self->fixedPar, par);
    ParStore_sanityCheck(self, __FILE__, __LINE__);
}

/// Add constrained parameter to ParStore.
void ParStore_addConstrainedPar(ParStore * self, const char *str,
                                const char *name, unsigned type) {
    Param *par = StrParMap_search(self->byname, name);
    if(par) {
        fprintf(stderr, "%s:%d: Duplicate definition of parameter \"%s\".\n",
                __FILE__, __LINE__, name);
        exit(EXIT_FAILURE);
    }

    int i = self->nPar++;
    if(i == MAXPAR) {
        fprintf(stderr, "%s:%s:%d: buffer overflow."
                " nPar=%d. MAXPAR=%d."
                " Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__, self->nPar, MAXPAR);
        exit(1);
    }

    par = self->vec + i;
    Param_init(par, name, 0.0, -DBL_MAX, DBL_MAX, type | CONSTRAINED);
    par->formula = strdup(str);
    CHECKMEM(par->formula);

    self->byname = StrParMap_insert(self->byname, par);
    self->byaddr = AddrParMap_insert(self->byaddr, par);
    self->constrainedPar = Param_push(self->constrainedPar, par);

    int status;
    par->constr = te_compile(str, self->te_pars, &status);
    if(par->constr == NULL) {
        fprintf(stderr, "%s:%d: parse error\n", __FILE__, __LINE__);
        fprintf(stderr, "  %s\n", str);
        fprintf(stderr, "  %*s^\nError near here\n", status - 1, "");
        exit(EXIT_FAILURE);
    }
    SET_CONSTR(par);
    self->te_pars = te_variable_push(self->te_pars, par->name, &par->value);
    ParStore_sanityCheck(self, __FILE__, __LINE__);
}

/// Get name of i'th free parameter
const char *ParStore_getNameFree(ParStore * self, int i) {
    assert(self);
    int j = 0;
    for(Param *par = self->freePar; par != NULL; par = par->next) {
        if(j == i)
            return par->name;
        ++j;
    }
    return NULL;
}

/// Return pointer associated with parameter name.
double *ParStore_findPtr(ParStore * self, unsigned * type,
                         const char *name) {
    Param *par = StrParMap_search(self->byname, name);
    if(par == NULL)
        return NULL;
    *type = par->type;
    return &par->value;
}

/// Return 1 if ptr is the address of a constrained parameter; 0
/// otherwise.  
int ParStore_isConstrained(const ParStore * self, const double *ptr) {
    assert(self);
    Param *par = AddrParMap_search(self->byaddr, ptr);
    if(par && (par->type & CONSTRAINED))
        return 1;
    return 0;
}

void ParStore_chkDependencies(ParStore * self, const double *valptr,
                              PtrSet * seen) {
    assert(self);

    if(self->byaddr == NULL)
        return;
    Param *par = AddrParMap_search (self->byaddr, valptr);
    if(par == NULL) {
        fprintf(stderr, "%s:%s:%d: can't find parameter with address %p\n",
                __FILE__, __func__, __LINE__, valptr);
        exit(EXIT_FAILURE);
    }

    // Nothing to be done unless par is constrained.
    if( !(par->type & CONSTRAINED) )
        return;

    // index of par in array
    int ipar = par - self->vec;
    assert(ipar >= 0);
    assert(ipar < MAXPAR);

    // Get list of pointers to parameters on which par depends.
    int len = 100;
    const double *dep[len];
    len = te_dependencies(par->constr, len, dep);

    // Check that each dependendant constrained parameter is in
    // "seen". This implies that dependencies will be set before par
    // is set.
    for(int j = 0; j < len; ++j) {
        // dpar is Param structure associated with value pointer dep[j]
        Param *dpar = AddrParMap_search(self->byaddr, dep[j]);
        if(dpar == NULL) {
            fprintf(stderr, "%s:%s:%d: can't find parameter with address %p\n",
                    __FILE__, __func__, __LINE__, dep[j]);
            exit(EXIT_FAILURE);
        }

        if( !(dpar->type & CONSTRAINED) ) {
            // dep[j] not constrained: no problem
            continue;
        }
        // Does par depend on itself?
        if(par == dpar) {
            fprintf(stderr, "%s:%d: Error: \"%s\" depends on itself\n",
                    __FILE__, __LINE__, par->name);
            exit(EXIT_FAILURE);
        }
        // If x and y are constrained parameters and y = f(x), then x
        // must come before y in the .lgo file.
        if(par < dpar) {
            fprintf(stderr, "%s:%d: Error: \"%s\" depends on"
                    " \"%s\" and must come later in\n"
                    " the .lgo file.\n",
                    __FILE__, __LINE__, par->name, dpar->name);
            exit(EXIT_FAILURE);
        }

        if(PtrSet_exists(seen, dep[j])) {
            // dep[i] is set before par: no problem. 
            continue;
        }
        // Pointer dep[j] is not in "seen", which means that it will
        // not be set before "par" is evaluated. Print an error
        // message and abort.
        fprintf(stderr, "%s:%d: Error: \"%s\" depends on \"%s\".\n",
                __FILE__, __LINE__, par->name, dpar->name);
        fprintf(stderr, "   When one constrained parameter depends on"
                " another, the dependent\n"
                "   parameter must appear in the tree, either earlier"
                " on the same\n"
                "   line of descent or on a different line of descent.\n");
        exit(EXIT_FAILURE);
    }
}

void ParStore_sanityCheck(ParStore * self, const char *file, int line) {
#  ifndef NDEBUG
    REQUIRE(self, file, line);
    REQUIRE(self->nPar >= 0, file, line);
    REQUIRE(self->nPar <= MAXPAR, file, line);

    // For each name: (1) make sure it's a legal name;
    // (2) get the pointer associated with that name.
    int i;
    Param *par, *par2;
    for(i = 0; i < self->nPar; ++i) {
        par = self->vec + i;
        Param_sanityCheck(par, file, line);
        par2 = StrParMap_search(self->byname, par->name);
        REQUIRE(par == par2, file, line);
        par2 = AddrParMap_search(self->byaddr, &par->value);
        REQUIRE(par == par2, file, line);
    }

    int npar = ParStore_nFree(self) + ParStore_nFixed(self) +
        ParStore_nConstrained(self);
    REQUIRE(npar == self->nPar, file, line);

    for(par = self->freePar; par; par = par->next)
        REQUIRE(par->type & FREE, file, line);

    for(par = self->fixedPar; par; par = par->next)
        REQUIRE(par->type & FIXED, file, line);

    for(par = self->constrainedPar; par; par = par->next)
        REQUIRE(par->type & CONSTRAINED, file, line);

#  endif
}

/// Return 1 if two ParStore objects are equal; 0 otherwise.
int ParStore_equals(ParStore * lhs, ParStore * rhs) {
    int i;
    if(lhs == rhs)
        return 1;
    if(lhs->nPar != rhs->nPar) {
        return 0;
    }
    if(ParStore_nFree(lhs) != ParStore_nFree(rhs)) {
        return 0;
    }
    if(ParStore_nConstrained(lhs) != ParStore_nConstrained(rhs)) {
        return 0;
    }
    for(i=0; i < lhs->nPar; ++i) {
        int c = Param_compare(lhs->vec+i, rhs->vec+i);
        if(c!=0)
            return 0;
    }
    if(ParStore_constrain(lhs)) {
        fprintf(stderr, "%s:%d: free parameters violate constraints\n",
                __FILE__, __LINE__);
    }
    if(ParStore_constrain(rhs)) {
        fprintf(stderr, "%s:%d: free parameters violate constraints\n",
                __FILE__, __LINE__);
    }
    return 1;
}

/// If ptr points to a constrained parameter, then set its value.
void ParStore_constrain_ptr(ParStore * self, double *ptr) {
    assert(self);
    if(self->byaddr == NULL)
        return;
    Param *par = AddrParMap_search(self->byaddr, ptr);
    assert(par);

    // If ptr isn't a constrained parameter, then return immediately
    if( !(par->type & CONSTRAINED) )
        return;

    // set value of constrained parameter
    SET_CONSTR(par);
}

/// First check to see that free parameters obey boundary constraints.
/// If not, then return 1. Otherwise, set values of all constrained
/// parameters and return 0.
int ParStore_constrain(ParStore * self) {
    Param *par;
    for(par = self->freePar; par; par = par->next) {
        if(par->value < par->low)
            return 1;
        if(par->value > par->high)
            return 1;
    }
    for(par = self->constrainedPar; par; par = par->next)
        SET_CONSTR(par);
    return 0;
}

/// Make sure Bounds object is sane.
void Bounds_sanityCheck(Bounds * self, const char *file, int line) {
#  ifndef NDEBUG
    REQUIRE(self, file, line);
    REQUIRE(self->lo_twoN >= 0.0, file, line);
    REQUIRE(self->lo_twoN < self->hi_twoN, file, line);
    REQUIRE(self->lo_t >= 0.0, file, line);
    REQUIRE(self->lo_t < self->hi_t, file, line);
#  endif
}

/// Return 1 if two Bounds objects are equal; 0 otherwise.
int Bounds_equals(const Bounds * lhs, const Bounds * rhs) {
    if(lhs == rhs)
        return 1;
    return lhs->lo_twoN == rhs->lo_twoN
        && lhs->hi_twoN == rhs->hi_twoN
        && lhs->lo_t == rhs->lo_t && lhs->hi_t == rhs->hi_t;
}
