#ifndef ARR_DIFFEV_H
#  define ARR_DIFFEV_H
#  define MAXPOP  500
#  define MAXDIM  35

#  include "typedefs.h"
#  include "state.h"
#  include <assert.h>
#  include <stdbool.h>
#  include <gsl/gsl_rng.h>

typedef struct TaskArg TaskArg;
typedef struct DiffEv DiffEv;
typedef struct DiffEvPar DiffEvPar;
typedef enum   DEStatus DEStatus;

// DE sets its status to Running on entry. If status still equals
// Running on return, something is wrong.
// On return, status should equal one of the following:
//  ReachedGoal: DE stopped because yspread <= ytol
//  FinishedIterations: DE stopped after completing all iterations
//  Interrupted: DE stopped on a signal (SIGINT or SIGTERM)
enum DEStatus {ReachedGoal, FinishedIterations, Interrupted, Running};

struct DiffEvPar {
    int         dim, ptsPerDim, refresh, strategy, nthreads, verbose;
    unsigned long seed;
    double      F, CR;
    void       *jobData;
    SimSched   *simSched;
    double      ytol;    // stop when yspread <= ytol
    void       *(*JobData_dup) (const void *);
    void        (*JobData_free) (void *);
    double      (*objfun) (int dim, double x[dim], void *, void *);

    // Set these equal to NULL unless you want each thread to maintain
    // state variables, which are passed to each job. This is useful,
    // for example, if you want to allocate a random number generator
    // that is used sequentially by all jobs executed by a single
    // thread.
    void       *threadData;
    void       *(*ThreadState_new) (void *);
    void        (*ThreadState_free) (void *);
    State      *state;
};

DEStatus diffev(int dim, double estimate[dim], double *loCost,
                double *yspread, DiffEvPar dep, gsl_rng * rng);
const char *diffEvStrategyLbl(int i);
void        sighandle(int signo);
void        handleSIGTERM(int signo);

#endif
