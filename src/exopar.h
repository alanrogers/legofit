#ifndef ARR_EXOPAR_H
#define ARR_EXOPAR_H

#include "typedefs.h"

double ExoPar_sample(ExoPar *self, double low, double high, gsl_rng *rng);
ExoParList *ExoParList_add(ExoParList *old, double *ptr, double m, double sd,
                           double sd, double low, double high);
void ExoParList_free(ExoParList *self);
ExoParTab *ExoParTab_new(ExoParList *list);
ExoPar const * const ExoParTab_find(ExoParTab *self, double *ptr)

#endif
