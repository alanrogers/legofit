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

#include "branchtab.h"
#include "cost.h"
#include "misc.h"
#include "network.h"
#include "patprob.h"
#include "simsched.h"
#include <gsl/gsl_rng.h>
#include <math.h>

/// Calculate cost.
/// @param[in] dim dimension of x
/// @param[in] x vector of parameter values.
/// @param jdata void pointer to a CostPar object, which contains
/// exogeneous parameters of the cost function.
/// @param tdata void pointer to a random number generator
/// @return cost
double costFun(int dim, double x[dim], void *jdata, void *tdata) {
    CostPar *cp = (CostPar *) jdata;
    gsl_rng *rng = (gsl_rng *) tdata;

    long nreps = SimSched_getSimReps(cp->simSched);
    DPRINTF(("%s:%d: nreps=%ld\n",__FILE__,__LINE__,nreps));

    if(Network_setParams(cp->network, dim, x))
        return HUGE_VAL;
    if(!Network_feasible(cp->network, 0))
        return HUGE_VAL;

    BranchTab  *prob = get_brlen(cp->network, nreps, cp->doSing, rng);
    BranchTab_normalize(prob);
#if COST==KL_COST
    double cost = BranchTab_KLdiverg(cp->obs, prob);
#elif COST==LNL_COST
    double cost = BranchTab_negLnL(cp->obs, prob);
#else
# error "Unknown cost method"
#endif

    BranchTab_free(prob);

    return cost;
}

/// Duplicate an object of class CostPar.
void * CostPar_dup(const void * arg) {
    assert(arg);
    const CostPar *old = (const CostPar *) arg;
    if(arg==NULL)
        return NULL;
    CostPar *new = memdup(old, sizeof(CostPar));
    CHECKMEM(new);
    new->obs = BranchTab_dup(old->obs);
    CHECKMEM(new->obs);
    new->network = Network_dup(old->network);
    CHECKMEM(new->network);
    new->simSched = old->simSched;
    CHECKMEM(new->simSched);
    return new;
}

/// CostPar destructor.
void CostPar_free(void *arg) {
    CostPar *self = (CostPar *) arg;
    BranchTab_free(self->obs);
    Network_free(self->network);
    if(self)
        free(self);
}
