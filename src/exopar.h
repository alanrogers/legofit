#ifndef ARR_EXOPAR_H
#define ARR_EXOPAR_H

#include "typedefs.h"
#include <gsl/gsl_rng.h>

void ExoParList_free(ExoParList *self);

ExoParTab *ExoParTab_new(void);
void ExoParTab_free(ExoParTab *self);
void ExoParTab_add(ExoParTab *self, double *ptr, double m, double sd,
                   double low, double high);
void ExoParTab_freeze(ExoParTab *self);
double ExoParTab_sample(ExoParTab *self, double *ptr,
                                double low, double high,
                                gsl_rng *rng);

#endif
