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
#include <stdbool.h>

struct ParStore {
    int         nFixed, nFree;  // number of parameters
    double      fixedVal[MAXPAR];   // parameter values
    double      freeVal[MAXPAR];    // parameter values
    double      loFree[MAXPAR]; // lower bounds
    double      hiFree[MAXPAR]; // upper bounds
    bool        time[MAXPAR];   // =true for time parameters
    char       *nameFixed[MAXPAR];  // Parameter names
    char       *nameFree[MAXPAR];   // Parameter names
    ParKeyVal  *head;           // linked list of name/ptr pairs
};

static int compareDblPtrs(const void *void_x, const void *void_y);
static int compareDbls(const void *void_x, const void *void_y);

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

void ParStore_print(ParStore *self, FILE *fp) {
    int i;
    fprintf(fp, "%5d fixed:\n", self->nFixed);
    for(i=0; i < self->nFixed; ++i)
        fprintf(fp, "   %8s = %lf\n", self->nameFixed[i], self->fixedVal[i]);
    ParStore_printFree(self, fp);
}

void ParStore_printFree(ParStore *self, FILE *fp) {
    int i;
    fprintf(fp, "%5d free :\n", self->nFree);
    for(i=0; i < self->nFree; ++i)
        fprintf(fp, "   %8s = %lf\n", self->nameFree[i], self->freeVal[i]);
}

/// Constructor
ParStore   *ParStore_new(void) {
    ParStore   *self = malloc(sizeof(ParStore));
    CHECKMEM(self);
    memset(self, 0, sizeof(ParStore));
    ParStore_sanityCheck(self, __FILE__, __LINE__);
    return self;
}

/// Duplicate a ParStore
ParStore   *ParStore_dup(const ParStore * old) {
    assert(old);
    ParStore   *new = memdup(old, sizeof(ParStore));
    new->head = NULL;

    int         i;
    for(i = 0; i < new->nFree; ++i) {
        new->nameFree[i] = strdup(old->nameFree[i]);
        new->head = ParKeyVal_add(new->head, new->nameFree[i],
                                  new->freeVal + i, true);
    }

    for(i = 0; i < new->nFixed; ++i) {
        new->nameFixed[i] = strdup(old->nameFixed[i]);
        new->head = ParKeyVal_add(new->head, new->nameFixed[i],
                                  new->fixedVal + i, false);
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

    ParKeyVal_free(self->head);
    free(self);
}

/// Add free parameter to ParStore.
void ParStore_addFreePar(ParStore * self, double value,
                         double lo, double hi, const char *name,
                         bool isTimeParameter) {
    int         i = self->nFree;

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
    self->time[i] = isTimeParameter;
    self->nameFree[i] = strdup(name);
    CHECKMEM(self->nameFree[i]);

    // Linked list associates pointer with parameter name.
    self->head = ParKeyVal_add(self->head, name, self->freeVal + i,
							   true);
}

/// Add fixed parameter to ParStore.
void ParStore_addFixedPar(ParStore * self, double value, const char *name) {
    int         i = self->nFixed;

    if(++self->nFixed >= MAXPAR)
        eprintf("%s:%s:%d: buffer overflow."
                " nFixed=%d. MAXPAR=%d."
                " Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__, self->nFixed, MAXPAR);

    self->fixedVal[i] = value;
    self->nameFixed[i] = strdup(name);
    CHECKMEM(self->nameFixed[i]);

    // Linked list associates pointer with parameter name.
    self->head = ParKeyVal_add(self->head, name, self->fixedVal + i,
		false);
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
double     *ParStore_findPtr(ParStore * self, bool *isfree,
							 const char *name) {
    return ParKeyVal_get(self->head, isfree, name);
}

void ParStore_sanityCheck(ParStore *self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self, file, line);
    REQUIRE(self->nFixed >= 0, file, line);
    REQUIRE(self->nFree >= 0, file, line);
    REQUIRE(self->nFixed < MAXPAR, file, line);
    REQUIRE(self->nFree < MAXPAR, file, line);

    // Make sure each name consists only of legal characters.
    int i;
    char *s;
    for(i=0; i < self->nFixed; ++i) {
        s = self->nameFixed[i];
        REQUIRE(NULL != s, file, line);
        REQUIRE(legalName(s), file, line);
    }
    for(i=0; i < self->nFree; ++i) {
        s = self->nameFree[i];
        REQUIRE(NULL != s, file, line);
        REQUIRE(legalName(s), file, line);
    }
    ParKeyVal_sanityCheck(self->head, file, line);
#endif    
}

int         ParStore_equals(const ParStore *lhs, const ParStore *rhs) {
    if(lhs == rhs)
        return 1;
    if(lhs->nFixed != rhs->nFixed)
        return 0;
    if(lhs->nFree != rhs->nFree)
        return 0;
    if(0 != memcmp(lhs->fixedVal, rhs->fixedVal,
                   lhs->nFixed*sizeof(lhs->fixedVal[0])))
        return 0;
    if(0 != memcmp(lhs->freeVal, rhs->freeVal,
                   lhs->nFree*sizeof(lhs->freeVal[0])))
        return 0;
    if(0 != memcmp(lhs->loFree, rhs->loFree,
                   lhs->nFree*sizeof(lhs->loFree[0])))
        return 0;
    if(0 != memcmp(lhs->hiFree, rhs->hiFree,
                   lhs->nFree*sizeof(lhs->hiFree[0])))
        return 0;
    int i;
    for(i=0; i < lhs->nFixed; ++i)
        if(0 != strcmp(lhs->nameFixed[i], rhs->nameFixed[i]))
            return 0;
    for(i=0; i < lhs->nFree; ++i)
        if(0 != strcmp(lhs->nameFree[i], rhs->nameFree[i]))
            return 0;
    return ParKeyVal_equals(lhs->head, rhs->head);
}

void Bounds_sanityCheck(Bounds * self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self, file, line);
    REQUIRE(self->lo_twoN >= 0.0, file, line);
    REQUIRE(self->lo_twoN < self->hi_twoN, file, line);
    REQUIRE(self->lo_t >= 0.0, file, line);
    REQUIRE(self->lo_t < self->hi_t, file, line);
#endif
}

int         Bounds_equals(const Bounds *lhs, const Bounds *rhs) {
    if(lhs == rhs)
        return 1;
    return lhs->lo_twoN == rhs->lo_twoN
        && lhs->hi_twoN == rhs->hi_twoN
        && lhs->lo_t == rhs->lo_t
        && lhs->hi_t == rhs->hi_t;
}

