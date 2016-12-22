#ifndef ARR_EXOPAR_H
#define ARR_EXOPAR_H

#include "typedefs.h"
#include <gsl/gsl_rng.h>

ExoPar     *ExoPar_new(void);
void        ExoPar_free(ExoPar * self);
void        ExoPar_add(ExoPar * self, double *ptr, double m, double sd);
void        ExoPar_freeze(ExoPar * self);
int         ExoPar_sample(ExoPar * self, double *ptr,
                          double low, double high, gsl_rng * rng);
ExoPar     *ExoPar_dup(const ExoPar * old);
void        ExoPar_shiftPtrs(ExoPar * self, size_t offset, int sign);

#endif
