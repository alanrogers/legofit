/**
 * @file point.c
 * @brief Class Point. Represents a point in parameter space, together with
 * its associated cost value.
 * @author Alan R. Rogers
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "point.h"
#include "misc.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int nPar = 0;
static size_t pointSize = 0;
pthread_mutex_t pointLock;

struct Point {
    double cost;
    double par[1]; // array of nPar parameter values, using struct hack
};

struct PointBuff {
    unsigned totpts;     // max number of points
    unsigned curpts;     // current number of entries in buffer
    unsigned buffbytes;  // number of bytes allocated in self->buf
    size_t curPtr;       // address where the next push will be stored
    size_t end;          // address just beyond end of PointBuff
    unsigned char buf[1]; // using struct hack
};

void Point_setNPar(int npar) {
    int status;

    status = pthread_mutex_lock(&pointLock);
    if(status)
        ERR(status, "lock");

    if(pointSize != 0) {
        fprintf(stderr, "%s:%d: %s can only be called once.\n",
                __FILE__,__LINE__,__func__);
        exit(EXIT_FAILURE);
    }
    nPar = npar;
    pointSize = sizeof(Point) + (npar-1)*sizeof(double);

    status = pthread_mutex_unlock(&pointLock);
    if(status)
        ERR(status, "unlock");
}

Point *Point_new(void) {
    if(pointSize == 0) {
        fprintf(stderr,"%s:%d: must call Point_setNPar before %s\n",
                __FILE__,__LINE__,__func__);
        exit(EXIT_FAILURE);
    }
    assert(nPar > 0);
    Point *self = malloc(pointSize);
    CHECKMEM(self);
    return self;
}

void Point_free(Point *self) {
    free(self);
}

void Point_set(Point *self, double cost, double *par) {
    assert(nPar > 0);
    self->cost = cost;
    memcpy(self->par, nPar, nPar*sizeof(par[0]));
}

double Point_get(Point *self, double *par) {
    assert(nPar > 0);
    memcpy(par, self->par, nPar*sizeof(par[0]));
    return cost;
}

PointBuff *PointBuff_new(unsigned totpts) {
    if(pointSize == 0) {
        fprintf(stderr,"%s:%d: must call Point_setNPar before %s\n",
                __FILE__,__LINE__,__func__);
        exit(EXIT_FAILURE);
    }

    self->buffbytes = totpts*pointsize;

    // Using struct hack.
    size_t size = sizeof(PointBuff) + self->buffbytes - 1;
    PointBuff *self = malloc(size);
    if(self==NULL)
        return NULL;

    self->end = ((size_t) self) + size;
    self->curPtr = &self->buff[0];
    self->curpts = 0;
    self->totpts = totpts;
    return self;
}

void PointBuff_free(Wraparound *self) {
    free(self);
}

// Return the number of items in the buffer.
unsigned PointBuff_size(const PointBuff *self) {
    return self->curpts;
}

void PointBuff_push(PointBuff *self, double cost, int n,
                    double param[n]) {
    Point *pt = (Point *) self->curPtr;
    assert(n == nPar);
    pt->cost = cost;
    memcpy(pt->par, param, nPar*sizeof(double));

    self->curPtr += pointSize;
    if(self->curPtr == self->end)
        self->curPtr = &self->buf[0];
    assert(self->curPtr < self->end);
    if(self->curpts != self->totpts)
        ++self->curpts;
}

// Put a vector of parameter values into param, return the corresponding
// value of cost, and remove the Point from which these values came
// from PointBuff.
double PointBuff_pop(PointBuff *self, int n, double param[n]) {
    if(self->curpts == 0) {
        fprintf(stderr,"%s:%d: can't pop an empty PointBufff\n",
                __FILE__, __LINE__);
        exit(1);
    }

    // Move to last filled position in buffer.
    self->curPtr = &self->buf[0] + pointSize*(self->curpts - 1);

    // Decrement number of points
    --self->curpts;

    // Return the point that has now been vacated
    Point *pt = (Point *) self->curPtr;
    memcpy(param, pt->par, nPar*sizeof(double));
    return pt->cost;
}
