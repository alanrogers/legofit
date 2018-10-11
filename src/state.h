#ifndef ARR_STATE_H
#define ARR_STATE_H

#include "typedefs.h"
#include <stdio.h>

NameList *NameList_append(NameList *self, const char *name);
void   NameList_free(NameList *self);
int    NameList_size(NameList *self);
void   NameList_print(NameList *self, FILE *fp);
int    State_npoints(State *self);
int    State_nparameters(State *self);
State *State_new(int npts, int npar);
void   State_free(State *self);
State *State_read(const char *fname, int npar, const char * const *name);
State *State_readList(NameList *list, int npts, int npar, const char * const *name);
int    State_print(State *self, FILE *fp);
int    State_setName(State *self, int ndx, const char *name);
void   State_setVector(State *self, int ndx, int dim, double x[dim]);
int    State_getVector(State *self, int ndx, int dim, double x[dim]);
void   State_setCost(State *self, int ndx, double cost);
double State_getCost(State *self, int ndx);

#endif
