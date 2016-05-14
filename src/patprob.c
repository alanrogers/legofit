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

typedef struct TaskArg TaskArg;

/** Data structure used by each thread */
struct TaskArg {
    unsigned    rng_seed;
    unsigned long nreps;
    GPTree     *gptree;

    // Returned value
    BranchTab  *branchtab;
};

TaskArg    *TaskArg_new(unsigned rng_seed, GPTree *gptree, unsigned nreps);
void        TaskArg_free(TaskArg * targ);
int         taskfun(void *varg);
static unsigned hashString(const char *ss);

/**
 * Construct a new TaskArg by copying a template, but then assign
 * a distinct random number seed.
 */
TaskArg    *TaskArg_new(unsigned rng_seed, GPTree *gptree, unsigned nreps) {
    TaskArg    *a = malloc(sizeof(TaskArg));
    checkmem(a, __FILE__, __LINE__);

    a->rng_seed = rng_seed;
    a->nreps = nreps;
    a->gptree = GPTree_dup(gptree);
    a->branchtab = BranchTab_new();

#if 0
    // Make sure all site patterns are represented in table.
    // This code iterates through all legal site patterns,
    // adding a branch of zero length for each pattern.
    tipId_t i;
    for(i=1; i < (1<<GPTree_nsamples(gptree))-1; ++i) {
        if(isPow2(i))
            continue;
        BranchTab_add(a->branchtab, i, 0.0);
    }
#endif

    return a;
}

/** TaskArg destructor */
void TaskArg_free(TaskArg * self) {
    BranchTab_free(self->branchtab);
    GPTree_free(self->gptree);
    free(self);
}

/** function run by each thread */
int taskfun(void *varg) {
    TaskArg    *targ = (TaskArg *) varg;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, targ->rng_seed);
    GPTree_simulate(targ->gptree, targ->branchtab, rng, targ->nreps);
    gsl_rng_free(rng);

    return 0;
}

/// Hash a character string
static unsigned hashString(const char *ss) {
    unsigned long hashval;
    int c;
    const unsigned char *s = (const unsigned char *) ss;
    // djb2
    hashval = 5381;
    while((c = *s++))
        hashval = ((hashval << 5) + hashval) +  c;

    return hashval;
}


/// Run simulations to estimate site pattern probabilities.
/// On return, pat[i] identifies the i'th pattern, and prob[i]
/// estimates its probability.  Function returns the number of
/// patterns detected--the length of pat and prob.
BranchTab *patprob(GPTree *gptree, int nThreads, long nreps, int pointNdx) {
    int j;
    unsigned long currtime = (unsigned long ) time(NULL);

    TaskArg    *taskarg[nThreads];
    unsigned    pid = (unsigned) getpid();

    /*
     * Generate a seed that is unique across points and threads.
     * First step creates a character string that concatenates
     * the decimal representation of time, process id, and index
     * of current point. This string is then hashed to a 32-bit
     * unsigned integer, which stored in a 64-bit unsigned
     * integer, lseed. Inside the loop, we add j to this seed to get
     * the seed for the j'th thread, and apply modulus to avoid
     * overflow.
     */
    long unsigned lseed;
    {
        char s[50]; // strlen(s) should not exceed 40
        snprintf(s, sizeof(s), "%lu%u%u", currtime, pid, pointNdx);
        assert(1+strlen(s) < sizeof(s));
        lseed = hashString(s);
    }

    // Divide repetitions among threads.
    long reps[nThreads];
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

    for(j = 0; j < nThreads; ++j) {
        unsigned seed = (unsigned) ((lseed + j) % UINT_MAX);
        taskarg[j] = TaskArg_new(seed, gptree, reps[j]);
    }

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

