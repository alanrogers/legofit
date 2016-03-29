#ifndef ARR_PARSTORE_H
#  define ARR_PARSTORE_H

#  include "typedefs.h"
#  include "misc.h"
#  include <assert.h>
#  define MAXPAR 50

struct Bounds {
    double      lo_twoN, hi_twoN, lo_t, hi_t;
};

ParStore   *ParStore_new(void);
void        ParStore_free(ParStore * self);
void        ParStore_addFreePar(ParStore * self, double value,
                                double lo, double hi, const char *name);
void        ParStore_addFixedPar(ParStore * self, double value,
                                 const char *name);
int         ParStore_nFixed(ParStore * self);
int         ParStore_nFree(ParStore * self);
double      ParStore_getFixed(ParStore * self, int i);
double      ParStore_getFree(ParStore * self, int i);
void        ParStore_setFree(ParStore * self, int i, double value);
double      ParStore_loFree(ParStore * self, int i);
double      ParStore_hiFree(ParStore * self, int i);
double     *ParStore_rawArray(ParStore * self);
double     *ParStore_findPtr(ParStore * self, const char *name);
ParStore   *ParStore_dup(ParStore *old);

#endif
