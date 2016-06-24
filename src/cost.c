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
    printf("%s:%s:%d\n",__FILE__,__func__,__LINE__); fflush(stdout);
    BranchTab_normalize(prob);
    printf("%s:%s:%d\n",__FILE__,__func__,__LINE__); fflush(stdout);
    double rval = BranchTab_KLdiverg(cp->obs, prob);
    printf("%s:%s:%d\n",__FILE__,__func__,__LINE__); fflush(stdout);
    return rval;
}
