#ifndef ARR_DIFFEV_H
#  define ARR_DIFFEV_H
#  define MAXPOP  500
#  define MAXDIM  35

#  include <assert.h>
#  include <stdbool.h>
#  include <gsl/gsl_rng.h>

typedef struct TaskArg TaskArg;
typedef struct DiffEv DiffEv;
typedef struct DiffEvPar DiffEvPar;

struct DiffEvPar {
    int         dim, ptsPerDim, genmax, refresh, strategy, nthreads, verbose;
    unsigned long seed;
    double      F, CR, deTol, costGoal;
    void       *jobData;
    double      (*objfun) (int dim, double x[dim], void *, void *);

    // Set these equal to NULL unless you want each thread to maintain
    // state variables, which are passed to each job. This is useful,
    // for example, if you want to allocate a random number generator
    // that is used sequentially by all jobs executed by a single
    // thread. 
    void       *threadData;
    void       *(*ThreadState_new) (void *);
    void        (*ThreadState_free) (void *);
    void       *initData;
    void        (*initialize)(int, void *, int, double *, gsl_rng *rng);
};

int         diffev(int dim, double estimate[dim], double *loCost,
                   double *yspread, DiffEvPar dep, gsl_rng * rng);
const char *diffEvStrategyLbl(int i);

#endif
