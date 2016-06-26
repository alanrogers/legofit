/**
 *        D I F F E R E N T I A L     E V O L U T I O N
 *
 * Program: de.c
 * Version: 3.6
 *
 * Authors: Dr. Rainer Storn
 *          c/o ICSI, 1947 Center Street, Suite 600
 *          Berkeley, CA 94707
 *          Tel.:   510-642-4274 (extension 192)
 *          Fax.:   510-643-7684
 *          E-mail: storn@icsi.berkeley.edu
 *          WWW: http://http.icsi.berkeley.edu/~storn/
 *          on leave from
 *          Siemens AG, ZFE T SN 2, Otto-Hahn Ring 6
 *          D-81739 Muenchen, Germany
 *          Tel:    636-40502
 *          Fax:    636-44577
 *          E-mail: rainer.storn@zfe.siemens.de
 *
 *          Kenneth Price
 *          836 Owl Circle
 *          Vacaville, CA 95687
 *          E-mail: kprice@solano.community.net
 *
 * This program implements some variants of Differential
 * Evolution (DE) as described in part in the techreport
 * tr-95-012.ps of ICSI. You can get this report either via
 * ftp.icsi.berkeley.edu/pub/techreports/1995/tr-95-012.ps.Z
 * or via WWW: http://http.icsi.berkeley.edu/~storn/litera.html
 * A more extended version of tr-95-012.ps is submitted for
 * publication in the Journal Evolutionary Computation.
 *
 * You may use this program for any purpose, give it to any
 * person or change it according to your needs as long as you
 * are referring to Rainer Storn and Ken Price as the originators of
 * the the DE idea. If you have questions concerning DE feel free to
 * contact us. We also will be happy to know about your experiences
 * with DE and your suggestions of improvement. 
 */

#include "diffev.h"
#include "misc.h"
#if 0
#include "jobqueue.h"
#endif
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#define DPRINTF_ON
#include "dprintf.h"
#ifdef DPRINTF_ON
extern pthread_mutex_t outputLock;
#endif

struct TaskArg {
    double      cost;
    int         dim;
    double     *v;
    double      (*objfun) (int dim, double x[dim], void *jdat, void *tdat);
    void       *jobData;
};

static inline void assignd(int dim, double a[], double b[]);
void        constrain(int n, double x[n], const double b[n],
                      const double loBnd[n], const double hiBnd[n]);
void        sample(int k, int rtn[k], int n, int array[n],
                   const gsl_rng * rng);
TaskArg    *TaskArg_new(int dim,
                        double (*objfun) (int xdim, double x[xdim],
                                          void *jdat, void *tdat),
                        void *jobData);
static inline void TaskArg_setArray(TaskArg * self, int dim, double v[dim]);
void        TaskArg_free(TaskArg * self);
void        TaskArg_print(TaskArg * self, FILE * fp);
int         taskfun(void *voidPtr, void *tdat);

static const char *stratLbl[] = // strategy-indicator
{ "", "DE/best/1/exp", "DE/rand/1/exp", "DE/rand-to-best/1/exp",
    "DE/best/2/exp", "DE/rand/2/exp", "DE/best/1/bin",
    "DE/rand/1/bin", "DE/rand-to-best/1/bin", "DE/best/2/bin",
    "DE/rand/2/bin"
};

/// Called by JobQueue
int taskfun(void *voidPtr, void *tdat) {
    TaskArg    *targ = (TaskArg *) voidPtr;
    targ->cost = targ->objfun(targ->dim, targ->v, targ->jobData, tdat);
    return 0;
}

/// If x violates boundary constraints, move it back inside.
/// This is done by shortening the vector (x - b), which represents
/// the direction of change. If x is changed, the revised x will reach
/// half way from b to the nearest boundary.
void constrain(int n, double x[n], const double b[n], const double loBnd[n],
               const double hiBnd[n]) {
    double      r, s = 1.0;
    int         i;

    for(i = 0; i < n; ++i) {
        assert(hiBnd[i] >= loBnd[i]);
        assert(b[i] >= loBnd[i]);
        assert(b[i] <= hiBnd[i]);

        if(x[i] > hiBnd[i]) {
            r = 0.5 * (hiBnd[i] - b[i]) / (x[i] - b[i]);
            s = fmin(r, s);
        } else if(x[i] < loBnd[i]) {
            r = 0.5 * (loBnd[i] - b[i]) / (x[i] - b[i]);
            s = fmin(r, s);
        }
        assert(r >= 0.0);
    }
    if(s < 1.0) {
        for(i = 0; i < n; ++i)
            x[i] = b[i] + s * (x[i] - b[i]);
    }
}

void TaskArg_print(TaskArg * self, FILE * fp) {
    int         i;
    fprintf(fp, "cost(");
    for(i = 0; i < self->dim; ++i) {
        fprintf(fp, "%lf", self->v[i]);
        if(i + 1 < self->dim)
            putc(',', fp);
        fprintf(fp, ")=%lf\n", self->cost);
    }
}

const char *diffEvStrategyLbl(int i) {
    unsigned    n = (unsigned) ((sizeof stratLbl) / (sizeof stratLbl[0]));
    if(i < 1 || i >= n) {
        fprintf(stderr, "%s:%s:%d: argument=%d; should be in [1, %u]\n",
                __FILE__, __func__, __LINE__, n, n-1u);
        exit(1);
    }
    return stratLbl[i];
}

TaskArg    *TaskArg_new(int dim,
                        double (*objfun) (int xdim, double x[xdim],
                                          void *jdat, void *tdat),
                        void *jobData) {
    TaskArg    *self = malloc(sizeof(TaskArg));
    CHECKMEM(self);

    self->cost = -1.0;
    self->dim = dim;
    self->jobData = jobData;
    self->objfun = objfun;
    self->v = malloc(dim * sizeof(self->v[0]));
    CHECKMEM(self->v);
    return self;
}

// Don't free self->jobData or self->objfun, which are not locally owned.
void TaskArg_free(TaskArg * self) {
    free(self->v);
    free(self);
}

#define PRARRAY(nn,xx) \
    do { \
      printf("%s:%s:%d:", __FILE__,__func__,__LINE__);    \
      int ii; \
      for(ii=0; ii < (nn); ++ii) \
          printf(" %lf", (xx)[ii]); \
      putchar('\n'); \
    }while(0);

static inline void TaskArg_setArray(TaskArg * self, int dim, double v[dim]) {
    assert(dim == self->dim);
    memcpy(self->v, v, dim * sizeof(v[0]));
    self->cost = -1.0;
}

/// Assigns dim-dimensional vector b to vector a.
static inline void assignd(int dim, double a[dim], double b[dim]) {
    memcpy(a, b, dim * sizeof(b[0]));
}

int diffev(int dim, double estimate[dim], double *loCost, double *yspread,
           DiffEvPar dep, gsl_rng * rng) {

    int         i, j, L, n;     // counting variables
    int         r[5];           // random indices
    int         imin;           // index to member with lowest energy
    int         gen;
    const int   genmax = dep.genmax;
    const int   refresh = dep.refresh;
    const int   strategy = dep.strategy;
    const double F = dep.F;
    const double CR = dep.CR;
    const double deTol = dep.deTol;
#if 0
    const int   nthreads = dep.nthreads;
#endif
    const int   verbose = dep.verbose;

    int         nPts = dep.dim * dep.ptsPerDim;
    int         status;

    double      loBound[dim];
    memcpy(loBound, dep.loBound, dim * sizeof(loBound[0]));

    double      hiBound[dim];
    memcpy(hiBound, dep.hiBound, dim * sizeof(hiBound[0]));

    int         ndx[nPts];
    for(i = 0; i < nPts; ++i)
        ndx[i] = i;

    double      c[nPts][dim], d[nPts][dim];

    *yspread = *loCost = strtod("NaN", NULL);

    double      tmp[dim], best[dim], bestit[dim];   // members
    double      cost[nPts];     // obj. funct. values
    double      cmin;           // help variables

#if 0
    JobQueue   *jq = JobQueue_new(nthreads, dep.threadData,
                                  dep.ThreadState_new,
                                  dep.ThreadState_free);
#endif

    TaskArg    *targ[nPts];
    printf("%s:%d initializing DE. nPts=%d\n",__FILE__,__LINE__,
           nPts); fflush(stdout);
    for(i = 0; i < nPts; ++i) {
        (*dep.randomize)(dep.randomizeData, dim, c[i], rng);
        targ[i] = TaskArg_new(dim, dep.objfun, dep.jobData);

        // calculate objective function values in parallel
        TaskArg_setArray(targ[i], dim, c[i]);
#if 0
        JobQueue_addJob(jq, taskfun, targ[i]);
#else
        taskfun(targ[i], NULL);
#endif
    }
#if 0
    JobQueue_waitOnJobs(jq);
#endif
    printf("%s:%d done initializing DE\n",__FILE__,__LINE__);
    fflush(stdout);
    cmin = HUGE_VAL;
    imin = INT_MAX;
    for(i = 0; i < nPts; i++) {
        cost[i] = targ[i]->cost;    // objective function value
		printf("%s:%d:cost[%d]=%lf\n",
			   __FILE__,__LINE__,i, cost[i]); fflush(stdout);
        if(cost[i] < cmin) {
            cmin = cost[i];
            imin = i;
        }
    }
    assert(imin < INT_MAX);
    assert(cmin < HUGE_VAL);
    assignd(dim, best, c[imin]);    // save best member ever
    assignd(dim, bestit, c[imin]);  // save best member of generation

    double      (*pold)[nPts][dim] = &c;    // old population (generation G)
    double      (*pnew)[nPts][dim] = &d;    // new population (generation G+1)
    double     *basevec = NULL;

    // Iteration loop
    for(gen = 1; gen <= genmax; ++gen) {
        printf("%s:%d: top of loop\n",__FILE__,__LINE__); fflush(stdout);
        imin = 0;
        for(i = 0; i < nPts; i++) { // Start of loop through ensemble

            // r gets indices of 5 random points.
            sample(5, r, nPts, ndx, rng);

            // copy i'th point into tmp
            assignd(dim, tmp, (*pold)[i]);

            // identify a random parameter to perturb first
            n = gsl_rng_uniform_int(rng, (unsigned long) dim);

            // Choice of strategy
            // We have tried to come up with a sensible
            // naming-convention: DE/x/y/z
            // DE :  stands for Differential Evolution
            // x  :  a string which denotes the vector to be perturbed
            // y  :  number of difference vectors taken for perturbation of x
            // z  :  crossover method (exp = exponential, bin = binomial)
            //
            // There are some simple rules which are worth following:
            // 1)  F is usually between 0.5 and 1 (in rare cases > 1)
            // 2)  CR is between 0 and 1 with 0., 0.3, 0.7 and
            //     1. being worth to be tried first=
            // 3)  To start off nPts = 10*dim is a reasonable
            //     choice. Increase nPts if misconvergence happens.
            // 4)  If you increase nPts, F usually has to be decreased
            // 5)  When the DE/best... schemes fail DE/rand... usually
            //     works and vice versa

            switch (strategy) {
            case 1:
                // DE/best/1/exp
                // Our oldest strategy but still not bad. However, we have
                // found several optimization problems where misconvergence
                // occurs.
                // strategy DE0 (not in our paper)
                L = 0;
                basevec = bestit;
                do {
                    tmp[n] =
                        bestit[n] + F * ((*pold)[r[1]][n] - (*pold)[r[2]][n]);
                    if(++n == dim)
                        n = 0;
                    ++L;
                } while((gsl_rng_uniform(rng) < CR) && (L < dim));
                break;
            case 2:
                // DE/rand/1/exp
                // This is one of my favourite strategies. It works
                // especially well when the "bestit[]"-schemes experience
                // misconvergence. Try e.g. F=0.7 and CR=0.5. as a first
                // guess. 
                // strategy DE1 in the techreport
                L = 0;
                basevec = (*pold)[r[0]];
                do {
                    tmp[n] =
                        (*pold)[r[0]][n] + F * ((*pold)[r[1]][n] -
                                                (*pold)[r[2]][n]);
                    if(++n == dim)
                        n = 0;
                    ++L;
                } while((gsl_rng_uniform(rng) < CR) && (L < dim));
                break;
            case 3:
                // DE/rand-to-best/1/exp
                // This strategy seems to be one of the best
                // strategies. Try F=0.85 and CR=1. If you get
                // misconvergence try to increase nPts. If this doesn't help
                // you should play around with all three control
                // variables. 
                // similiar to DE2 but generally better
                L = 0;
                basevec = (*pold)[i];
                do {
                    tmp[n] =
                        tmp[n] + F * (bestit[n] - tmp[n]) +
                        F * ((*pold)[r[0]][n] - (*pold)[r[1]][n]);
                    if(++n == dim)
                        n = 0;
                    ++L;
                } while((gsl_rng_uniform(rng) < CR) && (L < dim));
                break;
            case 4:
                // DE/best/2/exp is another powerful strategy worth trying
                L = 0;
                basevec = bestit;   // questionable choice
                do {
                    tmp[n] = bestit[n] +
                        F * ((*pold)[r[0]][n] + (*pold)[r[1]][n]
                             - (*pold)[r[2]][n] - (*pold)[r[3]][n]);
                    if(++n == dim)
                        n = 0;
                    ++L;
                } while((gsl_rng_uniform(rng) < CR) && (L < dim));
                break;
            case 5:
                // DE/rand/2/exp seems to be a robust optimizer for many
                // functions.
                L = 0;
                basevec = (*pold)[r[4]];    // questionable choice
                do {
                    tmp[n] = (*pold)[r[4]][n] +
                        F * ((*pold)[r[0]][n] + (*pold)[r[1]][n]
                             - (*pold)[r[2]][n] - (*pold)[r[3]][n]);
                    if(++n == dim)
                        n = 0;
                    ++L;
                } while((gsl_rng_uniform(rng) < CR) && (L < dim));
                break;
            case 6:
                // Remaining strategies have binomial crossover.
                // DE/best/1/bin
                basevec = bestit;
                for(L = 0; L < dim; ++L) {  // perform dim binomial trials 
                    if(L == 0 || gsl_rng_uniform(rng) < CR) {
                        tmp[n] = bestit[n] + F * ((*pold)[r[1]][n]
                                                  - (*pold)[r[2]][n]);
                    }
                    if(++n == dim)
                        n = 0;
                };
                break;
            case 7:
                // DE/rand/1/bin
                basevec = (*pold)[r[0]];
                for(L = 0; L < dim; ++L) {  // perform dim binomial trials
                    if(L == 0 || gsl_rng_uniform(rng) < CR) {
                        // change at least one parameter
                        tmp[n] =
                            (*pold)[r[0]][n] + F * ((*pold)[r[1]][n] -
                                                    (*pold)[r[2]][n]);
                    }
                    if(++n == dim)
                        n = 0;
                }
                break;
            case 8:
                // DE/rand-to-best/1/bin
                basevec = (*pold)[i];
                for(L = 0; L < dim; ++L) {  // perform dim binomial trials
                    if(L == 0 || gsl_rng_uniform(rng) < CR) {
                        // change at least one parameter
                        tmp[n] =
                            tmp[n] + F * (bestit[n] - tmp[n]) +
                            F * ((*pold)[r[0]][n] - (*pold)[r[1]][n]);
                    }
                    if(++n == dim)
                        n = 0;
                }
                break;
            case 9:
                // DE/best/2/bin
                basevec = bestit;   // questionable choice
                for(L = 0; L < dim; ++L) {  // perform dim binomial trials
                    if(L == 0 || gsl_rng_uniform(rng) < CR) {
                        // change at least one parameter
                        tmp[n] = bestit[n]
                            + F * ((*pold)[r[0]][n] + (*pold)[r[1]][n]
                                   - (*pold)[r[2]][n] - (*pold)[r[3]][n]);
                    }
                    if(++n == dim)
                        n = 0;
                }
                break;
            case 10:
                // DE/rand/2/bin
                basevec = (*pold)[r[4]];    // questionable choice
                for(L = 0; L < dim; ++L) {  // perform dim binomial trials
                    if(L == 0 || gsl_rng_uniform(rng) < CR) {
                        // change at least one parameter
                        tmp[n] = (*pold)[r[4]][n]
                            + F * ((*pold)[r[0]][n] + (*pold)[r[1]][n]
                                   - (*pold)[r[2]][n] - (*pold)[r[3]][n]);
                    }
                    if(++n == dim)
                        n = 0;
                }
                break;
            default:
                fprintf(stderr, "%s:%d: illegal strategy: %d\n",
                        __FILE__, __LINE__, strategy);
                exit(1);
            }
            constrain(dim, tmp, basevec, loBound, hiBound);

            // Trial mutation now in tmp[]. Calculate cost.
            TaskArg_setArray(targ[i], dim, tmp);
#if 0
            JobQueue_addJob(jq, taskfun, targ[i]);
#else
            printf("%s:%d: calling taskfun\n",__FILE__,__LINE__); fflush(stdout);
            taskfun(targ[i], NULL);
            printf("%s:%d: back from taskfun\n",__FILE__,__LINE__); fflush(stdout);
#endif
        }                       // End loop ensemble

#if 0
        JobQueue_waitOnJobs(jq);
#endif

        // 2nd pass through ensemble is the "reduce" portion
        // of "map-reduce". It generates a new generation, based
        // on the old generation and all the trials.
        double      cmax = -INFINITY;
        for(i = 0; i < nPts; ++i) {
            double      trial_cost = targ[i]->cost;
            if(trial_cost <= cost[i]) {
                // accept mutation
                cost[i] = trial_cost;
                assignd(dim, (*pnew)[i], targ[i]->v);
                if(trial_cost < cmin) { // Was this a new minimum? If so...
                    cmin = trial_cost;  // reset cmin to new low...
                    imin = i;
                    assignd(dim, best, targ[i]->v);
                }
            } else {
                // reject mutation: keep old value
                assignd(dim, (*pnew)[i], (*pold)[i]);
            }
            cmax = fmax(cmax, cost[i]);
        }

        // Best member of current generation
        assignd(dim, bestit, best);

        // swap population arrays. New generation becomes old one
        {
            double      (*pswap)[nPts][dim] = pold;
            pold = pnew;
            pnew = pswap;
        }

        // Difference between best and worst cost values
        *yspread = cmax - cmin;

        // Output part
        if(verbose && gen % refresh == 1) {
            // display after every refresh generations
            fprintf(stderr, "Generation=%d Objfun=%-15.10lg yspread=%lf\n",
                    gen, cmin, *yspread);
            fprintf(stderr, "   Best params:");
            for(j = 0; j < dim; j++) {
                fprintf(stderr, " %0.10lg", best[j]);
                if(j != dim - 1)
                    putc(',', stderr);
            }
            putc('\n', stderr);
        }
        printf("%s:%d bottom of loop\n",__FILE__,__LINE__); fflush(stdout);
        if(*yspread <= deTol)
            break;
    }
    // End iterations

#if 0
    JobQueue_noMoreJobs(jq);
#endif
    if(*yspread > deTol) {
        status = 1;
        if(verbose)
            fputs("No convergence\n", stdout);
    } else {
        status = 0;
        if(verbose)
            fputs("Converged\n", stdout);
    }

    // Return estimates
    *loCost = cmin;
    memcpy(estimate, best, dim * sizeof(estimate[0]));

    // Free memory
    for(i = 0; i < nPts; ++i)
        TaskArg_free(targ[i]);
#if 0
    JobQueue_free(jq);
#endif

    return status;
}

/// Choose at random k distinct integers from array of n, placing
/// result in "rtn". Function re-orders the entries of "array", but
/// leaves it otherwise unchanged. The function can therefore be
/// called repeatedly without re-initializing "array".
void sample(int k, int rtn[k], int n, int array[n], const gsl_rng * rng) {
    assert(k <= n);
    int         j;

    for(j = 0; j < k; ++j) {
        int         i = gsl_rng_uniform_int(rng, (unsigned long) n);
        rtn[j] = array[i];
        if(i != n - 1) {
            int         swap = array[i];
            array[i] = array[n - 1];
            array[n - 1] = swap;
        }
        --n;
    }
}
