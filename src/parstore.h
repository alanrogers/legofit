#ifndef ARR_PARSTORE_H
#  define ARR_PARSTORE_H

#  include "typedefs.h"
#  include "misc.h"
#  include <assert.h>
#  define MAXPAR 50

struct Bounds {	double lo_twoN, hi_twoN, lo_t, hi_t; };

ParStore   *ParStore_new(void);
void        ParStore_free(ParStore * self);
double     *ParStore_addPar(ParStore * self, int isfixed, double value,
							double lo, double hi, const char *name);
int ParStore_nFixed(ParStore * self);
int ParStore_nFixed(ParStore * self);
double ParStore_getFixed(ParStore * self, int i);
double ParStore_getFree(ParStore * self, int i);
void ParStore_setFree(ParStore * self, int i, double value);
double ParStore_loFixed(ParStore * self, int i);
double ParStore_hiFixed(ParStore * self, int i);
double ParStore_hiFree(ParStore * self, int i);
double * ParStore_getFreePtr(ParStore * self);
double *ParStore_getPtr(ParStore *self, const char *name);

#endif
