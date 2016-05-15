#ifndef PATPROB_INCLUDED
#  define PATPROB_INCLUDED

#include "typedefs.h"
BranchTab *patprob(int dim, double x[dim], const GPTree *gptree, int nThreads,
                   long nreps);
#endif
