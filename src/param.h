#ifndef LEGOFIT_PARAM_H
#define LEGOFIT_PARAM_H

#include "typedefs.h"
#include "tinyexpr.h"
#include <stdio.h>

#define NAMESIZE 40

struct Param {
    char *name;         // name of parameter; locally owned
    double value;
    double low, high;   // bounds;
    ParamType type;     // Free, Fixed, or Constrained
    char *formula;      // formula for constrained variable
    te_expr *constr;    // expression tree for constrained variable
    struct Param *next; // for a linked list of Param objects.
};

void   Param_init(Param *self, const char *name, double value,
                  double low, double high, ParamType type);
void Param_copy(Param *new, const Param *old);
Param *Param_push(Param *self, Param *new);
void   Param_freePtrs(Param *self);
void   Param_print(Param *self, FILE *fp);
int    Param_compare(const Param *lhs, const Param *rhs);

#endif
