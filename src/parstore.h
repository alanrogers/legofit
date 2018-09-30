#ifndef ARR_PARSTORE_H
#  define ARR_PARSTORE_H

#  include "typedefs.h"
#  include "misc.h"
#  include "gptree.h"
#  include <gsl/gsl_rng.h>
#  define MAXPAR 100

struct Bounds {
    double      lo_twoN, hi_twoN, lo_t, hi_t;
};

ParStore   *ParStore_new(void);
void        ParStore_free(ParStore * self);
void        ParStore_addFreePar(ParStore * self, double value, double lo,
                                double hi, const char *name, unsigned type);
void        ParStore_addFixedPar(ParStore * self, double value,
                                 const char *name, unsigned type);
void        ParStore_addConstrainedPar(ParStore * self, const char *str,
                                       const char *name, unsigned type);
int         ParStore_constrain(ParStore *self);
void        ParStore_constrain_ptr(ParStore *self, double *ptr);
int         ParStore_nFixed(ParStore * self);
int         ParStore_nFree(ParStore * self);
int         ParStore_nConstrained(ParStore * self);
const char *ParStore_getNameFree(ParStore * self, int i);
double     *ParStore_findPtr(ParStore * self, unsigned *type,
                             const char *name);
ParStore   *ParStore_dup(const ParStore * old);
void        ParStore_sanityCheck(ParStore * self, const char *file, int line);
void        ParStore_print(ParStore * self, FILE * fp);
void        ParStore_printFree(ParStore * self, FILE * fp);
void        ParStore_printConstrained(ParStore * self, FILE * fp);
int         ParStore_equals(ParStore * lhs, ParStore * rhs);
void        ParStore_setFreeParams(ParStore * self, int n, double x[n]);
void        ParStore_getFreeParams(ParStore * self, int n, double x[n]);
int         ParStore_isConstrained(const ParStore *self, const double *ptr);
void        ParStore_randomize(ParStore *self, GPTree *gpt, gsl_rng *rng);
void        Bounds_sanityCheck(Bounds * self, const char *file, int line);
int         Bounds_equals(const Bounds * lhs, const Bounds * rhs);

#endif
