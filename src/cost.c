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

/// KL divergence of observed relative to expected site pattern frequencies relative
///
/// @param[in] obs vector of observed site pattern frequencies
double costFun(int dim, double x[dim], void *jdata, void *notused) {
    const CostPar *cp = (CostPar *) jdata;

    // Kullback-Leibler divergence
    return KLdiverg(spdim, ofrq, efrq);
}
