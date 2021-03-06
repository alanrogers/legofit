#ifndef LEGOFIT_PARAM_H
#define LEGOFIT_PARAM_H

#include "typedefs.h"
#include "tinyexpr.h"
#include "strptrmap.h"
#include <stdio.h>
#include <gsl/gsl_randist.h>

#define NAMESIZE 40

struct Param {
    char *name;         // name of parameter
    double value;
    double low, high;   // bounds;
    unsigned type;      // (FREE | TWON) etc
    char *formula;      // formula for constrained variable
    te_expr *constr;    // expression tree for constrained variable
};

double Param_getTrialValue(const Param *self, gsl_rng *rng);
double Param_getValue(const Param *self);
int    Param_compare(const Param *lhs, const Param *rhs);
int    Param_isFree(const Param *self);
int    Param_setValue(Param *self, double value);
void   Param_compileConstraint(Param *self, StrPtrMap *te_pars);
void   Param_constrain(Param *par);
void   Param_copy(Param *to, const Param *from);
void   Param_freePtrs(Param *self);
void   Param_move(Param *to, Param *from);
Param *Param_new(const char *name, double value,
                 double low, double high,
                 unsigned type, const char *formula);
void   Param_sanityCheck(const Param *self, const char *file, int line);

#endif
