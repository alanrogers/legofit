#ifndef LEGOFIT_PARAM_H
#define LEGOFIT_PARAM_H

#include "typedefs.h"
#include <stdio.h>

#define FREE 1u
#define FIXED 2u
#define CONSTRAINED 4u
#define GAUSSIAN 8u
#define TWON 16u
#define TIME 32u
#define MIX 64u
#define NAMESIZE 40

struct Param {
    char *name;         // name of parameter; locally owned
    double value;
    double low, high;   // bounds;
    double mean, sd;    // parameters for gaussian variables
    unsigned type;     // combinations such as FIXED | TWON
};

void Param_init(Param *self, const char *name, double value,
                double low, double high,
                unsigned type);
void   Param_freePtrs(Param *self);
void   Param_print(Param *self, unsigned onlytype, FILE *fp);
void   Param_setGaussian(Param *self, double mean, double sd);

#endif
