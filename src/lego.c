/**
 * @file lego.c
 * @brief Simulate branch lengths
 *
 * @copyright Copyright (c) 2015, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "binary.h"
#include "branchtab.h"
#include "gptree.h"
#include "jobqueue.h"
#include "misc.h"
#include "parse.h"
#include "sampndx.h"
#include <assert.h>
#include <float.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

typedef struct TaskArg TaskArg;

/** Data structure used by each thread */
struct TaskArg {
    unsigned    rng_seed;
    unsigned long nreps;
    SampNdx     sndx;

    // Returned value
    BranchTab  *branchtab;
};

TaskArg    *TaskArg_new(const TaskArg * template, unsigned rng_seed);
void        TaskArg_free(TaskArg * targ);
int         taskfun(void *varg);

/**
 * Construct a new TaskArg by copying a template, but then assign
 * a distinct random number seed.
 */
TaskArg    *TaskArg_new(const TaskArg * template, unsigned rng_seed) {
    TaskArg    *a = malloc(sizeof(TaskArg));
    checkmem(a, __FILE__, __LINE__);

    memcpy(a, template, sizeof(TaskArg));
    a->rng_seed = rng_seed;

    return a;
}

/** TaskArg destructor */
void TaskArg_free(TaskArg * targ) {
    free(targ);
}

/** function run by each thread */
int taskfun(void *varg) {
    TaskArg    *targ = (TaskArg *) varg;

    int         i, bit[64];
    int         maxbits = sizeof(bit) / sizeof(bit[0]);
    unsigned long rep;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, targ->rng_seed);
    SampNdx     sndx;
    SampNdx_init(&sndx);
    BranchTab  *bt = BranchTab_new();

    PopNode    *rootPop = mktree(targ->fp, targ->ht, &sndx);

    for(rep = 0; rep < targ->nreps; ++rep) {
        PopNode_clear(rootPop);
        SampNdx_populateTree(&sndx);

        /* coalescent simulation */
        Gene       *root = PopNode_coalesce(rootPop, rng);
        assert(root);

        Gene_tabulate(root, targ->branchtab);
        Gene_free(root);
    }

    gsl_rng_free(rng);
    PopNode_free(rootPop);

    return 0;
}

int main(void) {

#if 0
    /* for production */
    int         nthreads = 28;  /* number of threads to launch */
    int         nTasks = 50;    /* total number of tasks */
    unsigned long nreps = 1000000;
#else
    /* for debugging */
    int         nthreads = 2;   /* number of threads to launch */
    int         nTasks = 5;     /* total number of tasks */
    unsigned long nreps = 100;
#endif

    TaskArg     targ = {
        .rng_seed = 0,
        .nreps = nreps,
    };

    TaskArg    *taskarg[nTasks];

    int         j;
    time_t      currtime = time(NULL);

    if(nthreads == 0)
        nthreads = getNumCores();

    if(nthreads > nTasks)
        nthreads = nTasks;

    printf("nreps       : %lu\n", nreps);
    printf("nthreads    : %d\n", nthreads);
    printf("nTasks      : %d\n", nTasks);

    for(j = 0; j < nTasks; ++j)
        taskarg[j] = TaskArg_new(&targ, (unsigned) currtime + j);

    JobQueue   *jq = JobQueue_new(nthreads);
    if(jq == NULL)
        eprintf("ERR@%s:%d: Bad return from JobQueue_new",
                __FILE__, __LINE__);
    for(j = 0; j < nTasks; ++j)
        JobQueue_addJob(jq, taskfun, taskarg[j]);
    JobQueue_waitOnJobs(jq);
    fprintf(stderr, "Back from threads\n");

    double      maxval = 0.0;
    printf("%8s %8s %8s %8s %8s %8s %8s %8s\n",
           "EQ", "oQ", "EIxy", "oIxy", "EInx", "oInx", "EIny", "oIny");
    for(j = 0; j < nTasks; ++j) {
        maxval = fmax(maxval, taskarg[j]->Qtheory);
        maxval = fmax(maxval, taskarg[j]->Qsim);
        printf("%8.4lf %8.4lf", taskarg[j]->Qtheory, taskarg[j]->Qsim);
        printf(" %8.4lf %8.4lf", taskarg[j]->Ixy,
               taskarg[j]->branchTab[1][0]);
        printf(" %8.4lf %8.4lf", taskarg[j]->Inx,
               taskarg[j]->branchTab[2][0]);
        printf(" %8.4lf %8.4lf", taskarg[j]->Iny,
               taskarg[j]->branchTab[2][1]);
        putchar('\n');
    }
    printf("\n%% max val: %0.4lf\n", maxval);

    for(j = 0; j < nTasks; ++j)
        TaskArg_free(taskarg[j]);

    return 0;
}
