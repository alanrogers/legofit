/**
   @file cost.c
   @page cost
   @brief Calculate cost function

   @copyright Copyright (c) 2016, Alan R. Rogers 
   <rogers@anthro.utah.edu>. This file is released under the Internet
   Systems Consortium License, which can be found in file "LICENSE".
*/

#undef DEBUG
#undef DPRINTF_ON

#include "dprintf.h"
#ifdef DPRINTF_ON
#include <pthread.h>
extern pthread_mutex_t outputLock;
#endif

#include "cost.h"
#include "gptree.h"
#include "branchtab.h"
#include "patprob.h"

/// @param[in] x vector of parameter values.
double costFun(int dim, double x[dim], void *jdata, void *notUsed) {
    const CostPar *cp = (CostPar *) jdata;

    BranchTab  *prob = patprob(dim, x, cp->gptree, cp->nThreads, cp->nreps);
    BranchTab_normalize(prob);
    double rval = BranchTab_KLdiverg(cp->obs, prob);
    return rval;
}
