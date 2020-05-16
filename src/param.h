#ifndef LEGOFIT_PARAM_H
#define LEGOFIT_PARAM_H

#include "typedefs.h"
#include "tinyexpr.h"
#include <stdio.h>
#include <gsl/gsl_randist.h>

#define NAMESIZE 40

struct Param {
    char *name;         // name of parameter; locally owned
    double value;
    double low, high;   // bounds;
    unsigned type;      // (FREE | TWON) etc
    char *formula;      // formula for constrained variable
    te_expr *constr;    // expression tree for constrained variable
    struct Param *next; // for a linked list of Param objects.
};

void   Param_init(Param *self, const char *name, double value,
                  double low, double high, unsigned type);
void   Param_copy(Param *new, const Param *old);
Param *Param_push(Param *self, Param *new);
void   Param_freePtrs(Param *self);
void   Param_print(Param *self, FILE *fp);
int    Param_compare(const Param *lhs, const Param *rhs);
void   Param_sanityCheck(const Param *self, const char *file, int line);
double Param_getTrialValue(Param *self, gsl_rng *rng);
int    Param_isFree(Param *self);
double Param_getValue(Param *self);
int    Param_setValue(Param *self, double value);

#endif
