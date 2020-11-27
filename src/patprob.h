#ifndef PATPROB_INCLUDED
#  define PATPROB_INCLUDED

#include "typedefs.h"
#include <gsl/gsl_rng.h>

BranchTab *get_brlen(const void *network, long nreps, int doSing,
                     double min_brlen, gsl_rng *rng);
#endif
