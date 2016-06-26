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
#include <math.h>

/// @param[in] x vector of parameter values.
double costFun(int dim, double x[dim], void *jdata, void *notUsed) {
	printf("%s:%s:%d: entry\n",__FILE__,__func__,__LINE__);
	fflush(stdout);
    const CostPar *cp = (CostPar *) jdata;

	if(!GPTree_feasible(cp->gptree)) {
		printf("%s:%d: tree isn't feasible\n",__FILE__,__LINE__);
		fflush(stdout);
		return HUGE_VAL;
	}

    BranchTab  *prob = patprob(dim, x, cp->gptree, cp->nThreads, cp->nreps);
    BranchTab_normalize(prob);
    double kl = BranchTab_KLdiverg(cp->obs, prob);
	printf("%s:%d: kl=%lf\n",__FILE__,__LINE__, kl);
	printf("%s:%s:%d: exit\n",__FILE__,__func__,__LINE__);
	fflush(stdout);
    return kl;
}
