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
#include "sampndx.h"
#include <stdlib.h>
#include <string.h>

typedef struct TaskArg TaskArg;

/** Data structure used by each thread */
struct TaskArg {
    const char *fname;
    unsigned    rng_seed;
    unsigned long nreps;
    SampNdx     sndx;
	Bounds      bnd;

    // Returned value
    BranchTab  *branchtab;
};

TaskArg    *TaskArg_new(const TaskArg * template, unsigned rng_seed,
                        unsigned nreps);
void        TaskArg_free(TaskArg * targ);
int         taskfun(void *varg);
char       *patLbl(size_t n, char buff[n], tipId_t tid, SampNdx * sndx);

/**
 * Construct a new TaskArg by copying a template, but then assign
 * a distinct random number seed.
 */
TaskArg    *TaskArg_new(const TaskArg * template, unsigned rng_seed,
                        unsigned nreps) {
    TaskArg    *a = malloc(sizeof(TaskArg));
    checkmem(a, __FILE__, __LINE__);

    memcpy(a, template, sizeof(TaskArg));
    a->rng_seed = rng_seed;
    a->nreps = nreps;
    a->branchtab = BranchTab_new();
    SampNdx_init(&(a->sndx));

    return a;
}

/** TaskArg destructor */
void TaskArg_free(TaskArg * self) {
    BranchTab_free(self->branchtab);
    free(self);
}

/** function run by each thread */
int taskfun(void *varg) {
    TaskArg    *targ = (TaskArg *) varg;

    unsigned long rep;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, targ->rng_seed);
    HashTab    *ht = HashTab_new();
	ParStore   *parstore = ParStore_new();  // parameters
    PopNode    *rootPop = NULL;
    {
        // Build population tree as specified in file targ->fname.
        // After this section, rootPop points to the ancestral
        // population, ht is a table that maps population names to
        // nodes in the population tree, and targ->sndx is an index of
        // samples. The call to HashTab_freeValues (at the end of this
        // function) deallocates all population nodes.
        FILE       *fp = fopen(targ->fname, "r");
        if(fp == NULL)
            eprintf("%s:%s:%d: can't open file %s.\n",
                    __FILE__, __func__, __LINE__, targ->fname);
        rootPop = mktree(fp, ht, &(targ->sndx), parstore, &(targ->bnd));
        fclose(fp);
    }

    for(rep = 0; rep < targ->nreps; ++rep) {
        PopNode_clear(rootPop); // remove old samples 
        SampNdx_populateTree(&(targ->sndx));    // add new samples

        // coalescent simulation generates gene genealogy within
        // population tree.
        Gene       *root = PopNode_coalesce(rootPop, rng);
        assert(root);

        // Traverse gene tree, accumulating branch lengths in bins
        // that correspond to site patterns.
        Gene_tabulate(root, targ->branchtab);

        // Free gene genealogy but not population tree.
        Gene_free(root);
    }

    gsl_rng_free(rng);
    HashTab_freeValues(ht);     // free all PopNode pointers
    HashTab_free(ht);
	ParStore_free(parstore);

    return 0;
}

/// Run simulations to estimate site pattern probabilities.
/// On return, pat[i] identifies the i'th pattern, and prob[i]
/// estimates its probability.  
unsigned patprob(unsigned maxpat,
                 tipId_t pat[maxpat],
                 double prob[maxpat],
                 int nTasks,
                 unsigned long nreps,
                 const char *fname,
                 Bounds *bnd,
                 unsigned rng_seed) {

    int i, j;

    // Divide repetitions among tasks.
    long        reps[nTasks];
    {
        ldiv_t      qr = ldiv((long) nreps, (long) nTasks);
        assert(qr.quot > 0);
        for(j = 0; j < nTasks; ++j)
            reps[j] = qr.quot;
        assert(qr.rem < nTasks);
        for(j=0; j < qr.rem; ++j)
            reps[j] += 1;
#ifndef NDEBUG
        // make sure the total number of repetitions is nreps.
        long        sumreps = 0;
        for(j = 0; j < nTasks; ++j) {
            assert(reps[j] > 0);
            sumreps += reps[j];
        }
        assert(sumreps = nreps);
#endif
    }

    TaskArg     targ = {
        .fname = fname,
        .rng_seed = 0,
        .nreps = 0,
		.bnd = bnd,
        .branchtab = NULL
    };

    TaskArg    *taskarg[nTasks];
    unsigned    pid = (unsigned) getpid();

    printf("# nreps       : %lu\n", nreps);
    printf("# nthreads    : %d\n", nTasks);
    printf("# input file  : %s\n", fname);

    for(j = 0; j < nTasks; ++j)
        taskarg[j] = TaskArg_new(&targ, currtime + pid + j, reps[j]);

    {
        JobQueue   *jq = JobQueue_new(nTasks);
        if(jq == NULL)
            eprintf("ERR@%s:%d: Bad return from JobQueue_new",
                    __FILE__, __LINE__);
        for(j = 0; j < nTasks; ++j)
            JobQueue_addJob(jq, taskfun, taskarg[j]);
        JobQueue_waitOnJobs(jq);
        JobQueue_free(jq);
    }
    fflush(stdout);

    // Add all branchtabs into branchtab[0]
    for(j = 1; j < nTasks; ++j)
        BranchTab_plusEquals(taskarg[0]->branchtab, taskarg[j]->branchtab);

    // Put site patterns and branch lengths into arrays.
    unsigned npat = BranchTab_size(taskarg[0]->branchtab);
    assert(npat <= maxpat);
    BranchTab_toArrays(taskarg[0]->branchtab, npat, pat, prob);

    {
        // Normalize so prob distribution sums to 1.
        double      sum = 0.0;
        for(j = 0; j < npat; ++j)
            sum += prob[j];
        for(j = 0; j < npat; ++j)
            prob[j] /= sum;
    }

    return npat;
}
