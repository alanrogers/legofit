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
#include "branchtab.h"
#include "parse.h"
#include "parstore.h"
#include "binary.h"
#include "gptree.h"
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
    GPTree     *gptree;

    // Returned value
    BranchTab  *branchtab;
};

static ThreadArg *ThreadArg_new(const GPTree *gptree, unsigned nreps,
                                int doSing);
static void ThreadArg_free(ThreadArg * targ);
static int tfunc(void *, void *);

/// function run by each thread
static int tfunc(void *varg, void *tdata) {
    ThreadArg    *arg = (ThreadArg *) varg;
    gsl_rng   *rng = (gsl_rng *) tdata;

    assert(GPTree_feasible(arg->gptree, 0));
    GPTree_patprobs(arg->gptree, arg->branchtab, rng, arg->nreps,
                    arg->doSing);

    return 0;
}

/// Construct a new ThreadArg by copying a template.
static ThreadArg *ThreadArg_new(const GPTree *gptree, unsigned nreps,
                            int doSing) {
    ThreadArg    *a = malloc(sizeof(ThreadArg));
    CHECKMEM(a);

    a->nreps = nreps;
    a->doSing = doSing;
    a->gptree = GPTree_dup(gptree);
    assert(GPTree_feasible(a->gptree, 0));
    a->branchtab = BranchTab_new();

    return a;
}

/// ThreadArg destructor
static void ThreadArg_free(ThreadArg * self) {
    BranchTab_free(self->branchtab);
    GPTree_free(self->gptree);
    free(self);
}

/// Estimate site pattern probabilities.  On return, pat[i] identifies
/// the i'th pattern, and prob[i] estimates its probability.  Function
/// returns a pointer to a newly-allocated object of type BranchTab,
/// which contains all the observed site patterns and their estimated
/// probabilities. 
BranchTab *patprob(const GPTree *gptree, long nreps,
                   int doSing, gsl_rng *rng) {

    ThreadArg *tharg;

    tharg = ThreadArg_new(gptree, nreps, doSing);
    tfunc(tharg, rng);

    BranchTab *rval = BranchTab_dup(tharg->branchtab);

    ThreadArg_free(tharg);

    return rval;
}
