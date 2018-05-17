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

struct Point {
    double cost;
    double par[1]; // array of nPar parameter values, using struct hack
};

Point *Point_new(void) {
    assert(nPar > 0);
    assert(pointSize > 0);
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
