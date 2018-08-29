#ifndef LEGOFIT_PARAM_H
#define LEGOFIT_PARAM_H

#include "typedefs.h"

struct Param {
    char *name;         // name of parameter; locally owned
    double *valptr;     // pointer to value; not locally owned
    double low, high;   // bounds;
    double mean, sd;    // parameters for gaussian variables
    ParamStatus status; // Free, Fixed, Gaussian, Constrained
};

Param *Param_new(const char *name, double *valptr, double low, double high,
                 ParamStatus status);
void   Param_free(Param *self);

#endif
