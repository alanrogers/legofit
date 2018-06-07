/**
 * @file pointbuff.c
 * @brief Class PointBuff. A buffer into which one can push points,
 * which consist of a vector of parameter values and an associated
 * cost value.
 * @author Alan R. Rogers
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "pointbuff.h"
#include "misc.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct Point {
    double cost;
    double par[1]; // array of nPar parameter values, using struct hack
};

struct PointBuff {
    unsigned nPar;
    size_t pointSize;
    unsigned totpts;     // max number of points
    unsigned curpts;     // current number of entries in buffer
    size_t curPtr;       // address where the next push will be stored
    size_t end;          // address just beyond end of PointBuff
    unsigned char buf[1]; // using struct hack
};

static inline void Point_set(Point *self, double cost, int n, double par[n]);
static inline double Point_get(Point *self, int n, double par[n]);

static inline void Point_set(Point *self, double cost, int n, double par[n]) {
    self->cost = cost;
    memcpy(self->par, par, n*sizeof(par[0]));
}

static inline double Point_get(Point *self, int n, double par[n]) {
    memcpy(par, self->par, n*sizeof(par[0]));
    return self->cost;
}

PointBuff *PointBuff_new(unsigned npar, unsigned totpts) {
    if(npar == 0) {
        fprintf(stderr,"%s:%d: can't allocate a PointBuff with 0 parameters\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    if(totpts == 0) {
        fprintf(stderr,"%s:%d: can't allocate a PointBuff with 0 points\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    size_t pointSize = sizeof(Point) + (npar-1)*sizeof(double);

    printf("sizeof(Point)=%zu pointSize = %zu\n",
           sizeof(Point), pointSize);

    // Using struct hack.
    size_t size = sizeof(PointBuff) + totpts*pointSize - 1;
    PointBuff *self = malloc(size);
    if(self==NULL)
        return NULL;

    self->nPar = npar;
    self->pointSize = pointSize;
    self->end = ((size_t) self) + size;
    self->curPtr = (size_t) self->buf;
    printf("Initialized curPtr=%zu end=%zu\n", self->curPtr, self->end);
    self->curpts = 0;
    self->totpts = totpts;
    return self;
}

void PointBuff_free(PointBuff *self) {
    free(self);
}

// Return the number of items in the buffer.
unsigned PointBuff_size(const PointBuff *self) {
    return self->curpts;
}

void PointBuff_push(PointBuff *self, double cost, int n,
                    double param[n]) {
    assert(n == self->nPar);
    Point *pt = (Point *) self->curPtr;
    Point_set(pt, cost, n, param);

    self->curPtr += self->pointSize;
    if(self->curPtr >= self->end) {
        printf("wrapping: diff=-%zu\n", self->curPtr - self->end);
        self->curPtr = (size_t) self->buf;
    }else{
        printf("diff=%zu\n", self->end - self->curPtr);
    }
    printf("curPtr=%zu end=%zu\n", self->curPtr, self->end);
    assert(self->curPtr < self->end);
    if(self->curpts != self->totpts)
        ++self->curpts;
}

// Put a vector of parameter values into param, return the corresponding
// value of cost, and remove the Point from which these values came
// from PointBuff.
double PointBuff_pop(PointBuff *self, int n, double param[n]) {
    assert(n == self->nPar);
    if(self->curpts == 0) {
        fprintf(stderr,"%s:%d: can't pop an empty PointBufff\n",
                __FILE__, __LINE__);
        exit(1);
    }

    // Move to last filled position in buffer.
    self->curPtr = (size_t) self->buf;
    self->curPtr += self->pointSize*(self->curpts - 1);

    // Decrement number of points
    --self->curpts;

    // Return the point that has now been vacated
    Point *pt = (Point *) self->curPtr;
    return Point_get(pt, self->nPar, param);
}
