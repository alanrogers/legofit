#ifndef ARR_PARSTORE_H
#  define ARR_PARSTORE_H

#  include "typedefs.h"
#  include "misc.h"
#  include <assert.h>
#  define MAXPAR 50

struct Bounds {	double lo_twoN, hi_twoN, lo_t, hi_t; };

ParStore   *ParStore_new(void);
void        ParStore_free(ParStore * self);
double     *ParStore_addPar(ParStore * self, double value, double lo, double hi);

int ParStore_nPar(ParStore * self);
double ParStore_get(ParStore * self, int i);
void ParStore_set(ParStore * self, int i, double value);
double ParStore_loBnd(ParStore * self, int i);
double ParStore_hiBnd(ParStore * self, int i);
double *ParStore_getPtr(ParStore * self);

#endif
