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
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "parstore.h"
#include "ptrqueue.h"
#include "ptrptrmap.h"
#include "strint.h"
#include "param.h"
#include "tinyexpr.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/errno.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <gsl/gsl_rng.h>

/// Abort with an error message about duplicate parameter definition
#define DUPLICATE_PAR(x) do{                                       \
        fprintf(stderr,"%s:%d: Duplicate parameter def: \"%s\"\n", \
                __FILE__,__LINE__, (x));                           \
        exit(EXIT_FAILURE);                                        \
    }while(0)

struct ParStore {

    int nPar, nFixed, nFree, nConstr; // nPar == nFixed+nFree+nConstr

    // An array of nPar Param objects. The 1st nFree entries describe
    // free parameters, the next nFixed describe fixed parameter
    // values, and the last nConstr entries describe constrained 
    // parameter values. I have them all in one array so that PopNode
    // objects can refer to them using a single index.
    Param *par;

    // These are pointers into the relevant portions of par.
    Param *free;    // equals par
    Param *fixed;   // equals par + nFree
    Param *constr;  // equals par + nFree + nFixed

    // Map parameter names to index numbers
    StrInt *byname;          // look up by name

    // Map the address of a parameter value to the address of the
    // enclosing Param object.
    PtrPtrMap *byaddr;      // look up by address

    te_variable *te_pars;    // for tinyexpr.c
};

static int  ParStore_constrain(ParStore *self);
static void ParStore_chkDependencies(ParStore * self);

ParStore *ParStore_new(PtrQueue *fixedQ, PtrQueue *freeQ, PtrQueue *constrQ) {
    ParStore *self = malloc(sizeof(ParStore));
    CHECKMEM(self);
    memset(self, 0, sizeof(ParStore));

    self->nFixed = PtrQueue_size(fixedQ);
    self->nFree = PtrQueue_size(freeQ);
    self->nConstr = PtrQueue_size(constrQ);
    self->nPar = self->nFixed + self->nFree + self->nConstr;

    self->par = malloc(self->nPar * sizeof(Param));
    CHECKMEM(self->par);

    self->free = self->par;
    self->fixed = self->free + self->nFree;
    self->constr = self->fixed + self->nFixed;

    // Param_move copies par into self->par.
    for(unsigned i=0; i < self->nFree; ++i) {
        Param *par = PtrQueue_pop(freeQ);
        Param_move(self->free+i, par);
        free(par);
    }

    for(unsigned i=0; i < self->nFixed; ++i) {
        Param *par = PtrQueue_pop(fixedQ);
        Param_move(self->fixed+i, par);
    }
    
    for(unsigned i=0; i < self->nConstr; ++i) {
        Param *par = PtrQueue_pop(constrQ);
        Param_move(self->constr+i, par);
    }

    self->byname = StrInt_new();
    self->byaddr = PtrPtrMap_new();
    self->te_pars = NULL;

    // Map parameter name to parameter index.
    // Map value address to Param address.
    // Push names and value pointers onto te_pars.
    for(unsigned i=0; i < self->nPar; ++i) {
        Param *par = self->par+i;

        errno = 0;
        StrInt_insert(self->byname, par->name, i);
        if(errno)
            DUPLICATE_PAR(par->name);

        int status = PtrPtrMap_insert(self->byaddr, &par->value, par);
        if(status) {
            fprintf(stderr,"%s:%d: duplicate parameter address\n",
                    __FILE__,__LINE__);
            exit(EXIT_FAILURE);
        }

        self->te_pars =
            te_variable_push(self->te_pars, par->name, &par->value);
    }

    // compile constraints
    for(unsigned i=0; i < self->nConstr; ++i)
        Param_compileConstraint(self->constr+i, self->te_pars);

    // Make sure no constrained parameter depends on a constrained
    // parameter that is defined later in the .lgo file.
    ParStore_chkDependencies(self);

    // set values of constrained variables
    if(ParStore_constrain(self)) {
        fprintf(stderr,"%s:%d:"
                " free parameters violate inequality constraints\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    
    ParStore_sanityCheck(self, __FILE__, __LINE__);
    return self;
}

/// Duplicate a ParStore
ParStore *ParStore_dup(const ParStore * old) {
    assert(old);
    ParStore *new = malloc(sizeof(ParStore));
    CHECKMEM(new);
    memset(new, 0, sizeof(ParStore));

    new->nFixed = old->nFixed;
    new->nFree = old->nFree;
    new->nConstr = old->nConstr;
    new->nPar = old->nPar;

    new->par = malloc(new->nPar * sizeof(Param));
    CHECKMEM(new->par);

    new->free = new->par;
    new->fixed = new->free + new->nFree;
    new->constr = new->fixed + new->nFixed;

    new->byname = StrInt_new();
    new->te_pars = NULL;

    for(int i = 0; i < new->nPar; ++i) {
        Param *par = new->par + i;

        const Param *opar = old->par + i;
        Param_copy(par, opar);  // doesn't copy constr field

        errno=0;
        StrInt_insert(new->byname, par->name, i);
        if(errno)
            DUPLICATE_PAR(par->name);

        new->te_pars =
            te_variable_push(new->te_pars, par->name, &par->value);
    }

    // compile constraints
    for(int i=0; i < new->nConstr; ++i)
        Param_compileConstraint(new->par+i, new->te_pars);

    ParStore_sanityCheck(new, __FILE__, __LINE__);
    return new;
}

/// Return the number of parameters
int ParStore_nPar(ParStore *self) {
    return self->nPar;
}

double ParStore_getVal(ParStore *self, int i) {
    assert(i >= 0);
    assert(i < self->nPar);
    return self->par[i].value;
}

Param *ParStore_getParam(ParStore *self, int i) {
    assert(i >= 0);
    assert(i < self->nPar);
    return self->par+i;
}


/// Return the number of free parameters
int ParStore_nFree(ParStore * self) {
    assert(self);
    return self->nFree;
}

/// Return the number of fixed parameters
int ParStore_nFixed(ParStore * self) {
    assert(self);
    return self->nFixed;
}

/// Return the number of constrained parameters
int ParStore_nConstrained(ParStore * self) {
    assert(self);
    return self->nConstr;
}

/// Set vector of free parameters, then update constrained parameters.
/// Return the value returned by ParStore_constrain.
int ParStore_setFreeParams(ParStore * self, int n, double x[n]) {
    assert(self);
    assert(n == self->nFree);
    for(int i=0; i < self->nFree; ++i)
        self->free[i].value = x[i];
    return ParStore_constrain(self);
}

/// Get vector of free parameters.
void ParStore_getFreeParams(ParStore * self, int n, double x[n]) {
    assert(self);
    assert(n == self->nFree);
    for(int i=0; i < self->nFree; ++i)
        x[i] = self->free[i].value;
}

/// Print a ParStore
void ParStore_print(ParStore * self, FILE * fp) {
    fprintf(fp, "Fixed:\n");
    for(int i=0; i < self->nFixed; ++i) 
        fprintf(fp, "   %8s = %lg\n", self->fixed[i].name,
                self->fixed[i].value);

    ParStore_printFree(self, fp);
    ParStore_printConstrained(self, fp);
}

/// Print free parameter values
void ParStore_printFree(ParStore * self, FILE * fp) {
    fprintf(fp, "Free:\n");
    for(int i=0; i < self->nFree; ++i)
        fprintf(fp, "   %8s = %lg\n", self->free[i].name, self->free[i].value);
}

/// Print constrained parameter values
void ParStore_printConstrained(ParStore * self, FILE * fp) {
    fprintf(fp, "Constrained:\n");
    for(int i=0; i < self->nConstr; ++i)
        fprintf(fp, "   %8s = %lg\n", self->constr[i].name,
                self->constr[i].value);
}

/// Destructor
void ParStore_free(ParStore * self) {
    assert(self);
    StrInt_free(self->byname);
    te_variable_free(self->te_pars);
    for(int i=0; i < self->nPar; ++i)
        Param_freePtrs(self->par + i);
    free(self->par);
    free(self);
}

/// Get name of i'th free parameter
const char *ParStore_getNameFree(ParStore * self, int i) {
    assert(self);
    assert(i < self->nFree);
    return self->free[i].name;
}

/// Return index of parameter name, or -1 if the name isn't there.
int ParStore_getIndex(ParStore * self, const char *name) {
    errno = 0;
    int i = StrInt_get(self->byname, name);
    if(errno)
        return -1;
    return i;
}

void ParStore_sanityCheck(ParStore * self, const char *file, int line) {
#  ifndef NDEBUG
    REQUIRE(self, file, line);
    REQUIRE(self->nPar >= 0, file, line);
    REQUIRE(self->nPar == self->nFree + self->nFixed + self->nConstr,
            file, line);
    
    // For each name: (1) make sure it's a legal name;
    // (2) get the pointer associated with that name.
    int i, j;
    Param *par;
    for(i = 0; i < self->nPar; ++i) {
        par = self->par + i;
        Param_sanityCheck(par, file, line);
        errno = 0;
        j = StrInt_get(self->byname, par->name);
        REQUIRE(errno==0, file, line);
        REQUIRE(j >= 0, file, line);
        REQUIRE(j == i, file, line);
        if(par->type & CONSTRAINED) {
            REQUIRE(par->formula, file, line);
        }else{
            REQUIRE(NULL == par->formula, file, line);
            REQUIRE(NULL == par->constr, file, line);
        }
    }

    for(i=0; i < self->nFixed; ++i)
        assert(FIXED & self->fixed[i].type);

    for(i=0; i < self->nFree; ++i)
        assert(FREE & self->free[i].type);

    for(i=0; i < self->nConstr; ++i)
        assert(CONSTRAINED & self->constr[i].type);
#  endif
}

/// Return 1 if two ParStore objects are equal; 0 otherwise.
int ParStore_equals(ParStore * lhs, ParStore * rhs) {
    if(lhs == rhs)
        return 1;
    if(lhs->nPar != rhs->nPar)
        return 0;
    if(ParStore_nFree(lhs) != ParStore_nFree(rhs))
        return 0;
    if(ParStore_nConstrained(lhs) != ParStore_nConstrained(rhs))
        return 0;
    for(int i=0; i < lhs->nPar; ++i) {
        if(0 != Param_compare(lhs->par+i, rhs->par+i) )
            return 0;
    }
    return 1;
}


/// First check to see that free parameters obey boundary constraints.
/// If not, then return 1. Otherwise, set values of all constrained
/// parameters and return 0.
static int ParStore_constrain(ParStore * self) {
    Param *par;
    for(int i=0; i < self->nFree; ++i) {
        par = self->free+i;
        if(par->value < par->low)
            return 1;
        if(par->value > par->high)
            return 1;
    }
    for(int i=0; i < self->nConstr; ++i)
        Param_constrain(self->constr + i);

    return 0;
}

static void ParStore_chkDependencies(ParStore * self) {
    assert(self);

    for(int i=0; i < self->nConstr; ++i ) {
        Param *par = self->constr + i;
        assert(par->type & CONSTRAINED);

        // Get list of pointers to parameters on which par depends.
        int len = 100;
        const double *dep[len];
        assert(par->constr != NULL);
        len = te_dependencies(par->constr, len, dep);

        // Check that each dependent constrained parameter comes before
        // "par" in the array of parameters. This implies that
        // dependencies will be set before par is set.
        for(int j = 0; j < len; ++j) {
            assert(self->byaddr);
            // dpar is Param structure associated with value pointer
            // dep[j]
            int status;
            Param *dpar = PtrPtrMap_get(self->byaddr, dep[j], &status);
            if(status) {
                fprintf(stderr, "%s:%s:%d:"
                        " can't find dependent parameter with address %p\n",
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
            // must be defined before y in the .lgo file. This implies
            // that par > dpar.
            if(par < dpar) {
                fprintf(stderr, "%s:%d: Error: \"%s\" depends on"
                        " \"%s\" and must be defined later in\n"
                        " the .lgo file.\n",
                        __FILE__, __LINE__, par->name, dpar->name);
                exit(EXIT_FAILURE);
            }
        }
    }
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
