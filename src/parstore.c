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
#include "parkeyval.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

struct ParStore {
    int         nFixed, nFree;      // number of parameters
    double      fixedVal[MAXPAR];   // parameter values
    double      freeVal[MAXPAR];    // parameter values
    double      loFree[MAXPAR];     // lower bounds
    double      hiFree[MAXPAR];     // upper bounds
    char       *nameFixed[MAXPAR];  // Parameter names
    char       *nameFree[MAXPAR];   // Parameter names
    ParKeyVal  *head;               // linked list of name/ptr pairs
};

/// Constructor
ParStore   *ParStore_new(void) {
    ParStore   *self = malloc(sizeof(ParStore));
    checkmem(self, __FILE__, __LINE__);
    memset(self, 0, sizeof(ParStore));
    return self;
}

/// Destructor
void ParStore_free(ParStore * self) {
    int         i;

    for(i = 0; i < self->nFixed; ++i)
        free(self->nameFixed[i]);

    for(i = 0; i < self->nFree; ++i)
        free(self->nameFree[i]);

    ParKeyVal_free(self->head);
    free(self);
}

/// Add free parameter to ParStore.
void ParStore_addFreePar(ParStore * self, double value,
                        double lo, double hi, const char *name) {
    int         i = self->nFree;

    if(++self->nFree >= MAXPAR)
        eprintf("%s:%s:%d: buffer overflow."
                " nFree=%d. MAXPAR=%d."
                " Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__,
                self->nFree, MAXPAR);

    if(value < lo || value > hi)
        eprintf("%s:%s:%d: value (%lf) not in range [%lf,%lf]\n",
                __FILE__, __func__, __LINE__,
                value, lo, hi);

    self->freeVal[i] = value;
    self->loFree[i] = lo;
    self->hiFree[i] = hi;
    self->nameFree[i] = strdup(name);
    CHECKMEM(self->nameFree[i]);

    // Linked list associates pointer with parameter name.
    self->head = ParKeyVal_add(self->head, name, self->freeVal + i);
}

/// Add fixed parameter to ParStore.
void ParStore_addFixedPar(ParStore * self, double value, const char *name) {
    int         i = self->nFixed;

    if(++self->nFixed >= MAXPAR)
        eprintf("%s:%s:%d: buffer overflow."
                " nFixed=%d. MAXPAR=%d."
                " Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__,
                self->nFixed, MAXPAR);

    self->fixedVal[i] = value;
    self->nameFixed[i] = strdup(name);
    CHECKMEM(self->nameFixed[i]);

    // Linked list associates pointer with parameter name.
    self->head = ParKeyVal_add(self->head, name, self->fixedVal + i);
}

/// Return the number of fixed parameters
int ParStore_nFixed(ParStore * self) {
    return self->nFixed;
}

/// Return the number of free parameters
int ParStore_nFree(ParStore * self) {
    return self->nFree;
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

/// Return pointer to array of free values
double     *ParStore_rawArray(ParStore * self) {
    return &self->freeVal[0];
}

/// Return pointer associated with parameter name.
double     *ParStore_findPtr(ParStore * self, const char *name) {
    return ParKeyVal_get(self->head, name);
}
