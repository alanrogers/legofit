/**
 * @file patprob.c
 * @author Alan R. Rogers
 * @brief Run simulations to estimate site pattern probabilities.
 *
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "patprob.h"
#include "misc.h"
#include "network.h"
#include "branchtab.h"
#include "parse.h"
#include "parstore.h"
#include "binary.h"
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>
#include <gsl/gsl_rng.h>

typedef struct ThreadArg ThreadArg;

/** Data structure used by each thread */
struct ThreadArg {
    unsigned long nreps;
    int         doSing; // nonzero => tabulate singletons
    void       *network;

    // Returned value
    BranchTab  *branchtab;
};

static ThreadArg *ThreadArg_new(const void *network, unsigned nreps,
                                int doSing);
static void ThreadArg_free(ThreadArg * targ);
static int tfunc(void *, void *);

/// function run by each thread
static int tfunc(void *varg, void *tdata) {
    ThreadArg    *arg = (ThreadArg *) varg;
    gsl_rng   *rng = (gsl_rng *) tdata;

    assert(Network_feasible(arg->network, 0));
    Network_brlen(arg->network, arg->branchtab, rng, arg->nreps,
                  arg->doSing);

    BranchTab_divideBy(arg->branchtab, (double) arg->nreps);

    return 0;
}

/// Construct a new ThreadArg by copying a template.
static ThreadArg *ThreadArg_new(const void *network, unsigned nreps,
                            int doSing) {
    ThreadArg    *a = malloc(sizeof(ThreadArg));
    CHECKMEM(a);

    a->nreps = nreps;
    a->doSing = doSing;
    a->network = Network_dup(network);
    assert(Network_feasible(a->network, 0));
    a->branchtab = BranchTab_new();

    return a;
}

/// ThreadArg destructor
static void ThreadArg_free(ThreadArg * self) {
    BranchTab_free(self->branchtab);
    Network_free(self->network);
    free(self);
}

/// Estimate branch lengths. Function returns a pointer to a
/// newly-allocated object of type BranchTab, which contains all the
/// observed site patterns and their mean branch lengths.
BranchTab *get_brlen(const void *network, long nreps,
                     int doSing, gsl_rng *rng) {

    ThreadArg *tharg;

    tharg = ThreadArg_new(network, nreps, doSing);
    tfunc(tharg, rng);

    BranchTab *rval = BranchTab_dup(tharg->branchtab);

    ThreadArg_free(tharg);

    return rval;
}
