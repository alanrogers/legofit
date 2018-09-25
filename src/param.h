#ifndef LEGOFIT_PARAM_H
#define LEGOFIT_PARAM_H

#include "typedefs.h"
#include <stdio.h>

#define NAMESIZE 40

struct Param {
    char *name;         // name of parameter; locally owned
    double value;
    double low, high;   // bounds;
    ParamType type;     // Free, Fixed, or Constrained
    char *formula;      // formula for constrained variable
    te_expr *constr;    // expression tree for constrained variable
};

void   Param_init(Param *self, const char *name, double value,
                  double low, double high,
                  ParamType type);
void   Param_freePtrs(Param *self);
void   Param_print(Param *self, ParamType onlytype, FILE *fp);

#endif
