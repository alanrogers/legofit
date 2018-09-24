#ifndef LEGOFIT_PARAM_H
#define LEGOFIT_PARAM_H

#include "typedefs.h"
#include <stdio.h>

#define FREE 1u
#define FIXED 2u
#define CONSTRAINED 4u
#define TWON 8u
#define TIME 16u
#define MIX 32u
#define GENERAL 64u
#define NAMESIZE 40

struct Param {
    char *name;         // name of parameter; locally owned
    double value;
    double low, high;   // bounds;
    unsigned type;      // combinations such as FIXED | TWON
    char *formula;      // formula for constrained variable
    te_expr *constr;    // expression tree for constrained variable
};

void   Param_init(Param *self, const char *name, double value,
                  double low, double high,
                  unsigned type);
void   Param_freePtrs(Param *self);
void   Param_print(Param *self, unsigned onlytype, FILE *fp);

#endif
