#ifndef ARR_PARSTORE_H
#  define ARR_PARSTORE_H

#  include "typedefs.h"
#  include "misc.h"
#  include <assert.h>
#  include <stdbool.h>
#  define MAXPAR 100

struct Bounds {
    double      lo_twoN, hi_twoN, lo_t, hi_t;
};

ParStore   *ParStore_new(void);
void        ParStore_free(ParStore * self);
void        ParStore_addFreePar(ParStore * self, double value,
                                double lo, double hi, const char *name);
void        ParStore_addGaussianPar(ParStore * self, double mean, double sd,
                                    const char *name);
void        ParStore_addFixedPar(ParStore * self, double value,
                                 const char *name);
void        ParStore_addConstrainedPar(ParStore * self, char *str,
                                       const char *name);
void        ParStore_constrain(ParStore *self);
void        ParStore_constrain_ptr(ParStore *self, double *ptr);
int         ParStore_nFixed(ParStore * self);
int         ParStore_nFree(ParStore * self);
int         ParStore_nGaussian(ParStore * self);
int         ParStore_nConstrained(ParStore * self);

double      ParStore_getFixed(ParStore * self, int i);
double      ParStore_getFree(ParStore * self, int i);
double      ParStore_getGaussian(ParStore * self, int i);

const char *ParStore_getNameFixed(ParStore * self, int i);
const char *ParStore_getNameFree(ParStore * self, int i);
const char *ParStore_getNameGaussian(ParStore * self, int i);

void        ParStore_setFree(ParStore * self, int i, double value);
double      ParStore_loFree(ParStore * self, int i);
double      ParStore_hiFree(ParStore * self, int i);
double     *ParStore_loBounds(ParStore * self);
double     *ParStore_upBounds(ParStore * self);
double     *ParStore_findPtr(ParStore * self, ParamStatus *pstat,
                             const char *name);
ParStore   *ParStore_dup(const ParStore * old);
void        ParStore_sanityCheck(ParStore * self, const char *file, int line);
void        ParStore_print(ParStore * self, FILE * fp);
void        ParStore_printFree(ParStore * self, FILE * fp);
void        ParStore_printConstrained(ParStore * self, FILE * fp);
int         ParStore_equals(const ParStore * lhs, const ParStore * rhs);
void        ParStore_setFreeParams(ParStore * self, int n, double x[n]);
void        ParStore_getFreeParams(ParStore * self, int n, double x[n]);
void        ParStore_sample(ParStore * self, double *ptr, double low,
                            double high, gsl_rng * rng);

void        Bounds_sanityCheck(Bounds * self, const char *file, int line);
int         Bounds_equals(const Bounds * lhs, const Bounds * rhs);

Constraint *Constraint_new(ParKeyVal * pkv, char *str);
Constraint *Constraint_free(Constraint * self);
double      Constraint_getValue(Constraint * self);
void        Constraint_prFormula(Constraint * self, FILE * fp);
Constraint *Constraint_dup(Constraint * old, ParKeyVal * pkv);
int         Constraint_equals(Constraint * lhs, Constraint * rhs);

#endif
