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
#  include <pthread.h>
extern pthread_mutex_t outputLock;
#endif

#include "cost.h"

/// @param[in] x vector of parameter values.
double costFun(int dim, double x[dim], void *jdata, void *tdata) {
    const CostPar *cp = (CostPar *) jdata;
	gsl_rng *rng = (gsl_rng *) tdata;

    GPTree_setParams(cp->gptree, dim, x);
    BranchTab *prob = patprob(cp->gptree, cp->nThreads, cp->nreps,
                              cp->pointNdx, rng);
    BranchTab_normalize(prob);
    return BranchTab_KLdiverg(prob, cp->obs);
}
