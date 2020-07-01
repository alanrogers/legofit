#ifndef ARR_PARSTORE_H
#  define ARR_PARSTORE_H

#  include "typedefs.h"
#  include "misc.h"
#  include <gsl/gsl_rng.h>
#  define MAXPAR 100

struct Bounds {
    double      lo_twoN, hi_twoN, lo_t, hi_t;
};

ParStore    ParStore_new(PtrQueue *fixedQ, PtrQueue *freeQ,
                         PtrQueue *constrQ);
void        ParStore_free(ParStore * self);
int         ParStore_nPar(ParStore *self);
double      ParStore_getVal(ParStore *self, int i);
int         ParStore_constrain(ParStore *self);
void        ParStore_constrain_ndx(ParStore *self, int i);
int         ParStore_nFixed(ParStore * self);
int         ParStore_nFree(ParStore * self);
int         ParStore_nConstrained(ParStore * self);
const char *ParStore_getNameFree(ParStore * self, int i);
int         ParStore_getIndex(ParStore * self, const char *name);
ParStore   *ParStore_dup(const ParStore * old);
void        ParStore_sanityCheck(ParStore * self, const char *file, int line);
void        ParStore_print(ParStore * self, FILE * fp);
void        ParStore_printFree(ParStore * self, FILE * fp);
void        ParStore_printConstrained(ParStore * self, FILE * fp);
int         ParStore_equals(ParStore * lhs, ParStore * rhs);
int         ParStore_setFreeParams(ParStore * self, int n, double x[n]);
void        ParStore_getFreeParams(ParStore * self, int n, double x[n]);
void        ParStore_chkDependencies(ParStore * self);
void        Bounds_sanityCheck(Bounds * self, const char *file, int line);
int         Bounds_equals(const Bounds * lhs, const Bounds * rhs);

#endif
