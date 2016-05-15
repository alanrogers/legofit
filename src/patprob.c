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

typedef struct TaskArg TaskArg;

/** Data structure used by each thread */
struct TaskArg {
    unsigned long nreps;
    GPTree     *gptree;

    // Returned value
    BranchTab  *branchtab;
};

TaskArg    *TaskArg_new(GPTree *gptree, unsigned nreps);
void        TaskArg_free(TaskArg * targ);
int         taskfun(void *varg);

/// function run by each thread
int taskfun(void *varg) {
    TaskArg    *targ = (TaskArg *) varg;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);

	// Lock seed, initialize random number generator, increment seed,
	// and unlock.
	pthread_mutex_lock(&seedLock);
    gsl_rng_set(rng, rngseed);
	rngseed = (rngseed == ULONG_MAX ? 0 : rngseed+1);
	pthread_mutex_unlock(&seedLock);

    GPTree_simulate(targ->gptree, targ->branchtab, rng, targ->nreps);
    gsl_rng_free(rng);

    return 0;
}

/// Construct a new TaskArg by copying a template.
TaskArg    *TaskArg_new(GPTree *gptree, unsigned nreps) {
    TaskArg    *a = malloc(sizeof(TaskArg));
    checkmem(a, __FILE__, __LINE__);

    a->nreps = nreps;
    a->gptree = GPTree_dup(gptree);
    a->branchtab = BranchTab_new();

    return a;
}

/// TaskArg destructor
void TaskArg_free(TaskArg * self) {
    BranchTab_free(self->branchtab);
    GPTree_free(self->gptree);
    free(self);
}

/// Run simulations to estimate site pattern probabilities.  On
/// return, pat[i] identifies the i'th pattern, and prob[i] estimates
/// its probability.  Function returns a pointer to a newly-allocated
/// object of type BranchTab, which contains all the observed site
/// patterns and their summed branch lengths.
BranchTab *patprob(GPTree *gptree, int nThreads, long nreps, int pointNdx) {
    int j;
    TaskArg    *taskarg[nThreads];
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

    for(j = 0; j < nThreads; ++j)
        taskarg[j] = TaskArg_new(gptree, reps[j]);

    {
        JobQueue   *jq = JobQueue_new(nThreads);
        if(jq == NULL)
            eprintf("s:%s:%d: Bad return from JobQueue_new",
                    __FILE__, __func__, __LINE__);
        for(j = 0; j < nThreads; ++j)
            JobQueue_addJob(jq, taskfun, taskarg[j]);
        JobQueue_waitOnJobs(jq);
        JobQueue_free(jq);
    }

    // Add all branchtabs into branchtab[0]
    for(j = 1; j < nThreads; ++j)
        BranchTab_plusEquals(taskarg[0]->branchtab, taskarg[j]->branchtab);

    BranchTab *rval = BranchTab_dup(taskarg[0]->branchtab);

    for(j = 0; j < nThreads; ++j)
        TaskArg_free(taskarg[j]);
        
    return rval;
}
