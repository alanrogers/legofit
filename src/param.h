#ifndef LEGOFIT_PARAM_H
#define LEGOFIT_PARAM_H

#include "typedefs.h"

#define FIXED 1
#define CONSTRAINED 2
#define GAUSSIAN 4
#define TWON 8
#define TIME 16
#define MIX 32
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
void   Param_free(Param *self);
void   Param_print(Param *self, unsigned onlytype; FILE *fp);
void   Param_setGaussianParams(Param *self, double mean, double sd);

#endif
