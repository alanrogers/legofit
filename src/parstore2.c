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
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "param.h"
#include "strparmap.h"
#include "addrparmap.h"
#include "parstore2.h"
#include <stdio.h>
#include <stdlib.h>


// Set value of i'th constraint. Abort with an error message if result is
// NaN.
#define SET_CONSTR(i) do {                                              \
        self->constrainedVal[(i)] = te_eval(self->constr[(i)]);         \
        if(isnan(self->constrainedVal[(i)])) {                          \
            fprintf(stderr,"%s:%d: constraint returned NaN\n",          \
                    __FILE__,__LINE__);                                 \
            fprintf(stderr,"formula: %s = %s\n",                        \
                    self->nameConstrained[(i)], self->formulas[(i)]);   \
            ParStore_printFree(self, stderr);                           \
            exit(EXIT_FAILURE);                                         \
        }                                                               \
    }while(0)

struct ParStore2 {
    int         nFixed, nFree, nGaussian, nConstrained; // num pars
    double      freeVal[MAXPAR];    // free parameter values
    double      fixedVal[MAXPAR];   // parameter values
    double      gaussianVal[MAXPAR]; // parameter values
    double      constrainedVal[MAXPAR]; // parameter values
    StrParMap  *byname; // look up parameters by name
    AddrParMap *byaddr; // look up parameters by address
    te_expr    *constr[MAXPAR];      // controls constrainedVal entries
    te_variable te_pars[MAXPAR];
    char       *formulas[MAXPAR];    // formulas of constrained vars

    // dependencies[i] records the names of the variables on which
    // the i'th constrained variable depends. Each constrained variable
    // must be defined in the .lgo after all its dependencies, and
    // constrained time variables cannot depend on constrained time
    // variables of descendants in the network of populations.
    StrInt     *dependencies[MAXPAR];
};

/// Return <0, 0, or >0, as x is <, ==, or > y.
int compareDblPtrs(const void *void_x, const void *void_y) {
    double * const * x = (double * const *) void_x;
    double * const * y = (double * const *) void_y;
    return **x - **y;
}

/// Return <0, 0, or >0, as x is <, ==, or > y.
int compareDbls(const void *void_x, const void *void_y) {
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
    Param *par;
    fprintf(fp, "%5d fixed:\n", self->nFixed);
    for(i=0; i < self->nFixed; ++i) {
        par = AddrParMap_search(self->byaddr, self->fixedVal+i);
        if(par==NULL)
            DIE("Can't find parameter address");
        fprintf(fp, "   %8s = %lg\n", par->name, self->fixedVal[i]);
    }
    fprintf(fp, "%5d Gaussian:\n", self->nGaussian);
    for(i=0; i < self->nGaussian; ++i) {
        par = AddrParMap_search(self->byaddr, self->gaussianVal+i);
        if(par==NULL)
            DIE("Can't find parameter address");
        fprintf(fp, "   %8s = Gaussian(%lg, %lg)\n", self->nameGaussian[i],
                par->mean, par->sd);
    }
    ParStore_printFree(self, fp);
    ParStore_printConstrained(self, fp);
}

