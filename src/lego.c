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
    unsigned rng_seed;
    unsigned long nreps;
    SampNdx sndx;

    // Returned value
    BranchTab *branchtab;
};

TaskArg *TaskArg_new(const TaskArg *template, unsigned rng_seed);
void TaskArg_free(TaskArg *targ);
int taskfun(void *varg);

/**
 * Construct a new TaskArg by copying a template, but then assign
 * a distinct random number seed.
 */
TaskArg *TaskArg_new(const TaskArg *template, unsigned rng_seed) {
    TaskArg *a = malloc(sizeof(TaskArg));
    checkmem(a, __FILE__, __LINE__);

    memcpy(a, template, sizeof(TaskArg));
    a->rng_seed = rng_seed;

    return a;
}

/** TaskArg destructor */
void TaskArg_free(TaskArg *targ) {
    free(targ);
}

/** function run by each thread */
int taskfun(void *varg) {
    TaskArg *targ = (TaskArg *) varg;

    int         i, bit[64];
    int         maxbits = sizeof(bit) / sizeof(bit[0]);
    unsigned long rep;
    PopNode    *X, *Y0, *N0, *Y1, *Y2, *N1, *D, *XY, *ND, *URHUMAN;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, targ->rng_seed);

    /* perturb parameters */
    if(targ->mD > 0.0)
        targ->mD = perturb_ratio(targ->mD, rng);

    targ->alpha = perturb_interval(targ->alpha, 0.5*targ->alpha, targ->delta, rng);
    targ->delta = perturb_interval(targ->delta, targ->alpha, targ->epsilon, rng);
    targ->epsilon = perturb_interval(targ->epsilon, targ->delta, targ->zeta, rng);
    targ->zeta = perturb_interval(targ->zeta, targ->epsilon, targ->kappa, rng);
    targ->kappa = perturb_interval(targ->kappa, targ->zeta, targ->lambda, rng);
    targ->lambda = perturb_interval(targ->lambda, targ->kappa, 2.0*targ->lambda, rng);
    assert(targ->alpha < targ->delta);
    assert(targ->delta < targ->epsilon);
    assert(targ->epsilon < targ->zeta);
    assert(targ->zeta < targ->kappa);
    assert(targ->kappa < targ->lambda);

    targ->K_X = perturb_ratio(targ->K_X, rng);
    targ->K_Y = perturb_ratio(targ->K_Y, rng);
    targ->K_D = perturb_ratio(targ->K_D, rng);
    targ->K_XY = perturb_ratio(targ->K_XY, rng);
    targ->K_ND = perturb_ratio(targ->K_ND, rng);
    targ->K_N = perturb_ratio(targ->K_N, rng);

    X = PopNode_new(targ->K_X, 0.0, targ->zeta);
    Y0 = PopNode_new(targ->K_Y, 0.0, targ->alpha);    /* before D mixture */
    N1 = PopNode_new(targ->K_N, targ->epsilon, targ->kappa); /* after fossil */
    N0 = PopNode_new(targ->K_N, targ->delta, targ->epsilon);  /* before fossil */
    Y1 = PopNode_new(targ->K_Y, targ->alpha, targ->delta); /* btw D and N mix */
    Y2 = PopNode_new(targ->K_Y, targ->delta, targ->zeta);  /* after N mixture */
    D = PopNode_new(targ->K_D, targ->alpha, targ->kappa);
    XY = PopNode_new(targ->K_XY, targ->zeta, targ->lambda);
    ND = PopNode_new(targ->K_ND, targ->kappa, targ->lambda);
    URHUMAN = PopNode_new(1.0, targ->lambda, DBL_MAX);
    
    /* Y0 is a mixture: mD*D + (1-mD)*Y1 */
    PopNode_mix(Y0, targ->mD, D, Y1);

    /* Y1 is a mixture: mN*N0 + (1-mN)*Y2 */
    PopNode_mix(Y1, targ->mN, N0, Y2);

    /* N1 is the immediate ancestor of N0 */
    PopNode_endToEnd(N0, N1);

    PopNode_join(XY, X, Y2);
    PopNode_join(ND, N1, D);
    PopNode_join(URHUMAN, XY, ND);

    for(rep = 0; rep < targ->nreps; ++rep) {
        PopNode_clear(URHUMAN);
        PopNode_newGene(X, 0);
        PopNode_newGene(Y0, 1);
        PopNode_newGene(N1, 2);

        /* coalescent simulation */
        Gene       *root = PopNode_coalesce(URHUMAN, rng);
        assert(root);

        double      branch;
        tipId_t     tipId;

        Gene_tabulate(root, targ->branchtab);
        Gene_free(root);
    }

#ifndef NDEBUG
    {
        unsigned long sumcount = 0;
        int j;

        for(i=0; i < 3; ++i)
            for(j=0; j<i; ++j)
            sumcount += targ->count[i][j];
        assert(sumcount == targ->nreps);
    }
#endif

    gsl_rng_free(rng);
    PopNode_free(X);
    PopNode_free(Y0);
    PopNode_free(N0);
    PopNode_free(Y1);
    PopNode_free(Y2);
    PopNode_free(N1);
    PopNode_free(D);
    PopNode_free(XY);
    PopNode_free(ND);
    PopNode_free(URHUMAN);

    return 0;
}

int main(void) {

    double      twoN0 = 2e4;        /* haploid size of ancestral human pop */
    double      gen = 29.0;         /* generation time */
    double      s = twoN0 * gen;    /* years per time unit */

#if 0
    /* for production */
    int         nthreads = 28;   /* number of threads to launch */
    int         nTasks = 50;     /* total number of tasks */
    unsigned long   nreps = 1000000;
#else
    /* for debugging */
    int         nthreads = 2;   /* number of threads to launch */
    int         nTasks = 5;     /* total number of tasks */
    unsigned long   nreps = 10000;
#endif

    TaskArg targ = {
        .mN = 0.02,          /* N->EV gene flow */
        .mD = 0.03,          /* D->V gene flow */

        /* Time backwards from the present, units of twoN0 gen */
        .alpha = 25e3/s,     /* Denisovan admixture */
        .delta = 55e3/s,     /* Neanderthal admixture */
        .epsilon = 65e3/s,      /* age of older Neanderthal fossil */
        .zeta = 110e3/s,     /* Y and EV split */
        .kappa = 427e3/s,    /* N and D split */
        .lambda = 658e3/s,   /* archaics and moderns split */

        /* population sizes relative to N0 */
        .K_X = 1.0,
        .K_Y = 1.0,
        .K_D = 1.0,
        .K_XY = 1.0,
        .K_ND = 1.0,
        .K_N = 1.0,

        .rng_seed = 0,
        .nreps = nreps,

        /* Returned values */
        .count = {{0}},
        .branchTab = {{0.0}},
        .Inx = 0.0,
        .Iny = 0.0,
        .Ixy = 0.0,
        .Qsim = 0.0,
        .Qtheory = 0.0
    };

    TaskArg *taskarg[nTasks];

#if 1
    /* Populations differ in size */
    targ.K_X = 2.0;
    targ.K_XY = 2.0;
    targ.K_ND = 0.2;
    targ.K_N = 0.1;
    targ.K_D = 0.1;
#elif 0
    /* Force coalescent events into ancestral population */
    targ.K_X = targ.K_Y = targ.K_D = targ.K_XY = targ.K_ND =
        targ.K_N = strtod("Inf", 0);
#endif

    int         j;
    time_t      currtime = time(NULL);

    if(nthreads == 0)
        nthreads = getNumCores();

    if(nthreads > nTasks)
        nthreads = nTasks;

    printf("Q\n");
    printf("2N0         : %0.0lf\n", twoN0);
    printf("generation  : %0.1lf y\n", gen);
    printf("t/2N0       : delta=%lf alpha=%lf epsilon=%lf\n",
           targ.delta, targ.alpha, targ.epsilon);
    printf("            : kappa=%lf lambda=%lf zeta=%lf\n",
           targ.kappa, targ.lambda, targ.zeta);
    printf("Rel pop size: X=%lf Y=%lf N=%lf D=%lf\n",
           targ.K_X, targ.K_Y, targ.K_N, targ.K_D);
    printf("            : XY=%lf ND=%lf URHUMAN=%lf\n",
           targ.K_XY, targ.K_ND, 1.0);
    printf("admixture->Y: mN=%lf mD=%lf\n",
           targ.mN, targ.mD);
    printf("nreps       : %lu\n", nreps);
    printf("nthreads    : %d\n", nthreads);
    printf("nTasks      : %d\n", nTasks);

    for(j = 0; j < nTasks; ++j)
        taskarg[j] = TaskArg_new(&targ, (unsigned) currtime+j );

    JobQueue   *jq = JobQueue_new(nthreads);
    if(jq == NULL)
        eprintf("ERR@%s:%d: Bad return from JobQueue_new",
                __FILE__, __LINE__);
    for(j = 0; j < nTasks; ++j)
        JobQueue_addJob(jq, taskfun, taskarg[j]);
    JobQueue_waitOnJobs(jq);
    fprintf(stderr, "Back from threads\n");

    double maxval = 0.0;
    printf("%8s %8s %8s %8s %8s %8s %8s %8s\n",
           "EQ", "oQ", "EIxy", "oIxy", "EInx", "oInx", "EIny", "oIny");
    for(j=0; j<nTasks; ++j) {
        maxval = fmax(maxval, taskarg[j]->Qtheory);
        maxval = fmax(maxval, taskarg[j]->Qsim);
        printf("%8.4lf %8.4lf", taskarg[j]->Qtheory, taskarg[j]->Qsim);
        printf(" %8.4lf %8.4lf", taskarg[j]->Ixy, taskarg[j]->branchTab[1][0]);
        printf(" %8.4lf %8.4lf", taskarg[j]->Inx, taskarg[j]->branchTab[2][0]);
        printf(" %8.4lf %8.4lf", taskarg[j]->Iny, taskarg[j]->branchTab[2][1]);
        putchar('\n');
    }
    printf("\n%% max val: %0.4lf\n", maxval);

    double x[nTasks], y[nTasks];

    /* plot observed versus expected Q */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Qtheory;
        y[j] = taskarg[j]->Qsim;
    }
    pictex(x, y, nTasks, "EQ", "oQ", "$Q$", "figQsim.tex");

    /* plot observed versus expected Ixy */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Ixy;
        y[j] = taskarg[j]->branchTab[1][0];
    }
    pictex(x, y, nTasks, "EIxy", "oIxy", "$I_{xy}/ML$",
           "figQIxy.tex");

    /* plot observed versus expected Inx */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Inx;
        y[j] = taskarg[j]->branchTab[2][0];
    }
    pictex(x, y, nTasks, "EInx", "oInx", "$I_{nx}/ML$",
           "figQInx.tex");

    /* plot observed versus expected Iny */
    for(j=0; j<nTasks; ++j) {
        x[j] = taskarg[j]->Iny;
        y[j] = taskarg[j]->branchTab[2][1];
    }
    pictex(x, y, nTasks, "EIny", "oIny", "$I_{ny}/ML$",
           "figQIny.tex");

    for(j=0; j<nTasks; ++j)
        TaskArg_free(taskarg[j]);
    
    return 0;
}
