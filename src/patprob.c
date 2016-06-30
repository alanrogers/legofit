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
#include "jobqueue.h"
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

#undef DPRINTF_ON
#include "dprintf.h"
#ifdef DPRINTF_ON
extern pthread_mutex_t outputLock;
#endif

typedef struct SimArg SimArg;

/** Data structure used by each thread */
struct SimArg {
    unsigned long nreps;
    GPTree     *gptree;

    // Returned value
    BranchTab  *branchtab;
};

SimArg     *SimArg_new(const GPTree *gptree, unsigned nreps);
void        SimArg_free(SimArg * targ);
int         simfun(void *, void *);

/// function run by each thread
int simfun(void *varg, void *notUsed) {
    SimArg    *targ = (SimArg *) varg;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);

	// Lock seed, initialize random number generator, increment seed,
	// and unlock.
	pthread_mutex_lock(&seedLock);
    gsl_rng_set(rng, rngseed);
	rngseed = (rngseed == ULONG_MAX ? 0 : rngseed+1);
	pthread_mutex_unlock(&seedLock);

	assert(GPTree_feasible(targ->gptree));
    GPTree_simulate(targ->gptree, targ->branchtab, rng, targ->nreps);
    gsl_rng_free(rng);

    return 0;
}

/// Construct a new SimArg by copying a template.
SimArg    *SimArg_new(const GPTree *gptree, unsigned nreps) {
    SimArg    *a = malloc(sizeof(SimArg));
    checkmem(a, __FILE__, __LINE__);

    a->nreps = nreps;
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
BranchTab *patprob(const GPTree *gptree, int nThreads, long nreps) {

	assert(GPTree_feasible(gptree));
	
    int j;
    SimArg    *simarg[nThreads];
    long reps[nThreads];

    // Divide repetitions among threads.
    {
        ldiv_t      qr = ldiv(nreps, (long) nThreads);
        assert(qr.quot > 0);
        for(j = 0; j < nThreads; ++j)
            reps[j] = qr.quot;
        assert(qr.rem < nThreads);
        for(j=0; j < qr.rem; ++j)
            reps[j] += 1;
#ifndef NDEBUG
        // make sure the total number of repetitions is nreps.
        long        sumreps = 0;
        for(j = 0; j < nThreads; ++j) {
            assert(reps[j] > 0);
            sumreps += reps[j];
        }
        assert(sumreps == nreps);
#endif
    }

	assert(GPTree_feasible(gptree));
    for(j = 0; j < nThreads; ++j)
        simarg[j] = SimArg_new(gptree, reps[j]);

    {
        JobQueue   *jq = JobQueue_new(nThreads, NULL, NULL, NULL);
        if(jq == NULL)
            eprintf("s:%s:%d: Bad return from JobQueue_new",
                    __FILE__, __func__, __LINE__);
        for(j = 0; j < nThreads; ++j)
            JobQueue_addJob(jq, simfun, simarg[j]);
        JobQueue_waitOnJobs(jq);
        JobQueue_free(jq);
    }

    // Add all branchtabs into branchtab[0]
    for(j = 1; j < nThreads; ++j)
        BranchTab_plusEquals(simarg[0]->branchtab, simarg[j]->branchtab);

    BranchTab *rval = BranchTab_dup(simarg[0]->branchtab);

    for(j = 0; j < nThreads; ++j)
        SimArg_free(simarg[j]);
        
    return rval;
}
