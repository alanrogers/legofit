/**
 * @file patprob.c
 * @brief Run simulations to estimate site pattern probabilities.
 */

#include "patprob.h"
#include "misc.h"
#include "branchtab.h"
#include "hashtab.h"
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

typedef struct SimArg SimArg;

/** Data structure used by each thread */
struct SimArg {
    unsigned long nreps;
    int         doSing; // nonzero => tabulate singletons
    GPTree     *gptree;

    // Returned value
    BranchTab  *branchtab;
};

SimArg     *SimArg_new(const GPTree *gptree, unsigned nreps, int doSing);
void        SimArg_free(SimArg * targ);
int         simfun(void *, void *);

/// function run by each thread
int simfun(void *varg, void *tdata) {
    SimArg    *arg = (SimArg *) varg;
    gsl_rng   *rng = (gsl_rng *) tdata;

	assert(GPTree_feasible(arg->gptree));
    GPTree_simulate(arg->gptree, arg->branchtab, rng, arg->nreps,
                    arg->doSing);

    return 0;
}

/// Construct a new SimArg by copying a template.
SimArg    *SimArg_new(const GPTree *gptree, unsigned nreps, int doSing) {
    SimArg    *a = malloc(sizeof(SimArg));
    checkmem(a, __FILE__, __LINE__);

    a->nreps = nreps;
    a->doSing = doSing;
    a->gptree = GPTree_dup(gptree);
	assert(GPTree_feasible(a->gptree));
    a->branchtab = BranchTab_new();

    return a;
}

/// SimArg destructor
void SimArg_free(SimArg * self) {
    BranchTab_free(self->branchtab);
    GPTree_free(self->gptree);
    free(self);
}

/// Run simulations to estimate site pattern probabilities.  On
/// return, pat[i] identifies the i'th pattern, and prob[i] estimates
/// its probability.  Function returns a pointer to a newly-allocated
/// object of type BranchTab, which contains all the observed site
/// patterns and their summed branch lengths.
BranchTab *patprob(const GPTree *gptree, int nThreads, long nreps,
                   int doSing, gsl_rng *rng) {

    SimArg    *simarg;

    simarg = SimArg_new(gptree, nreps, doSing);
    simfun(simarg, rng);

    BranchTab *rval = BranchTab_dup(simarg->branchtab);

    SimArg_free(simarg);

    return rval;
}
