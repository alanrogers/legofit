#ifndef PATPROB_INCLUDED
#  define PATPROB_INCLUDED

#include "typedefs.h"
BranchTab *patprob(const GPTree *gptree, int nThreads, long nreps,
                   int doSing);
#endif
