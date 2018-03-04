#ifndef ARR_STATE_H
#define ARR_STATE_H

#include "typedefs.h"
#include <stdio.h>

int    State_npoints(State *self);
int    State_nparameters(State *self);
State *State_new(int npts, int npar);
void   State_free(State *self);
State *State_read(FILE *fp);
int    State_print(State *self, FILE *fp);
void   State_setVector(State *self, int ndx, int dim, double x[dim]);
int    State_getVector(State *self, int ndx, int dim, double x[dim]);
void   State_setCost(State *self, int ndx, double cost);
double State_getCost(State *self, int ndx);

#endif
