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
#include <pthread.h>

pthread_mutex_t seedLock = PTHREAD_MUTEX_INITIALIZER;
unsigned long rngseed;

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
int simfun(void *varg, void *notUsed) {
    SimArg    *arg = (SimArg *) varg;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);

	// Lock seed, initialize random number generator, increment seed,
	// and unlock.
	pthread_mutex_lock(&seedLock);
    gsl_rng_set(rng, rngseed);
	rngseed = (rngseed == ULONG_MAX ? 0 : rngseed+1);
	pthread_mutex_unlock(&seedLock);

	assert(GPTree_feasible(arg->gptree));
    GPTree_simulate(arg->gptree, arg->branchtab, rng, arg->nreps,
                    arg->doSing);
    gsl_rng_free(rng);

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
                   int doSing) {

	assert(GPTree_feasible(gptree));

    SimArg    *simarg;

	assert(GPTree_feasible(gptree));
    simarg = SimArg_new(gptree, nreps, doSing);
    simfun(simarg, NULL);

    BranchTab *rval = BranchTab_dup(simarg->branchtab);

    return rval;
}
