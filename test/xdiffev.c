#include "diffev.h"
#include "simsched.h"
#include "state.h"
#include "pointbuff.h"
#include "misc.h"
#include <getopt.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <memory.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void        usage(void);
double      objFunc(int dim, double x[dim], void *jdat, void *tdat);
void        initState(int ndx, void *void_p, int n, double x[n],
                         gsl_rng *rng);
void prOpt(const char *opt, const char *description);

#define RUGGED

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

/// Local minima wherever all entries of v are integers. Global
/// minimum at v=0.
///
/// For constrained optimization, have objFunc return HUGE_VAL when
/// constraints are violated.
double objFunc(int dim, double x[dim], void *jdat /* NOTUSED */ ,
               void *tdat /* NOTUSED */ ) {
    int         i;
    double      cost, sx = 0.0;
#ifdef RUGGED
    // for multiple peaks
    double      sf = 0.0;
    double      f;              // fractional part of x
#endif

    for(i = 0; i < dim; ++i) {
        double      xi = x[i];
        sx += fabs(xi);         // summed absolute devs from zero
#ifdef RUGGED
        // for multiple peaks
        f = xi - floor(xi + 0.5);   // fractional part of x[i]
        sf += fabs(f);          // summed fractional dev
#endif
    }
    cost = 1.0 + sx;
#ifdef RUGGED
    // for multiple peaks
    cost *= 1.0 + sf;
#endif
    return cost;
}

/// Initialize vector x. If ndx==0, simply copy the parameter vector
/// from the argument. Otherwise, generate a random vector.
///
/// For constrained optimization, make all initial vectors obey constraints.
void initState(int ndx, void *void_p, int n, double x[n], gsl_rng *rng){
	assert(void_p != NULL);
    double *v = (double *) void_p; // pointer to n-vector
    if(ndx == 0)
		memcpy(x, v, n*sizeof(x[0]));
    else {
		int i;
		for(i=0; i < n; ++i)
			x[i] = gsl_ran_gaussian(rng, 5.0);
    }
}

#undef RUGGED

/// Describe an option. For use in "usage" functions.
/// @param[in] opt Name of option.
/// @param[in] description Description of option.
void prOpt(const char *opt, const char *description) {
    fprintf(stderr, "   %s\n      %s\n", opt, description);
    return;
}

/// Print usage message and exit.
void usage(void) {
    fprintf(stderr, "usage: diffev [options]\n");
    fprintf(stderr, "   where options may include:\n");

    /* misc */
    prOpt("-s <method> or --strategy <method>", "strategy");
    prOpt("-g <x> or --genmax <x>", "max generations");
    prOpt("-r <x> or --refresh <x>", "refresh interval");
    prOpt("-n or --nParam", " number of parameters");
    prOpt("-p <x> or --ptsPerDim <x>", "points per dimension");
    prOpt("-F <x> or --F <x>", "DE weight factor");
    prOpt("-c <x> or --crossOver <x>", "crossover probability");
    prOpt("-v or --verbose", "more output");
    prOpt("-h or --help", "print this message");
    exit(1);
}

int main(int argc, char *argv[]) {
    int         optndx;
    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
        {"strategy", required_argument, 0, 's'},
        {"genmax", required_argument, 0, 'g'},
        {"refresh", required_argument, 0, 'r'},
        {"nParam", required_argument, 0, 'n'},
        {"ptsPerDim", required_argument, 0, 'p'},
        {"F", required_argument, 0, 'F'},
        {"crossOver", required_argument, 0, 'c'},
        {"tol", required_argument, 0, 'T'},
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, 0, 'v'},
        {NULL, 0, NULL, 0}
    };

    int         verbose = 0;
    int         strategy = 1;   // which flavor of differential evolution
    int         dim = 2;        // Dimension of parameter vector
    int         ptsPerDim = 10; // points per dimension
    int         nPts = 0;       // number of population members
    int         genmax = 1000;
    int         refresh = 10;   // refresh rate
    double      F = 0.9;        // scale perturbations
    double      CR = 0.8;       // crossover prob
    double      ytol = 1e-4;    // convergence criterion

    int         i;
    int         nthreads = 2;
    time_t      currtime = time(NULL);
    unsigned    baseSeed = currtime % UINT_MAX;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    assert(rng);
    gsl_rng_set(rng, baseSeed);
    long simreps = 1000;
    SimSched    *simSched = SimSched_new();

    // command line arguments
    for(;;) {
        i = getopt_long(argc, argv, "T:s:g:r:n:p:F:c:hv", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            usage();
            break;
        case 's':
            strategy = strtol(optarg, NULL, 10);
            break;
        case 'T':
            ytol = strtod(optarg, 0);
            break;
        case 'g':
            genmax = strtol(optarg, NULL, 10);
            break;
        case 'r':
            refresh = strtol(optarg, NULL, 10);
            break;
        case 'n':
            dim = strtol(optarg, NULL, 10);
            break;
        case 'p':
            ptsPerDim = strtol(optarg, NULL, 10);
            break;
        case 'F':
            F = strtod(optarg, NULL);
            break;
        case 'c':
            CR = strtod(optarg, NULL);
            break;
        case 'v':
            verbose = 1;
            break;
        case 'h':
        default:
            usage();
        }
    }

    SimSched_append(simSched, 200, 1000);
    SimSched_append(simSched, 100, 10000);
    SimSched_append(simSched, 1000, simreps);

    // There should be no non-flag arguments
    if(argc > optind)
        usage();


	double      initVec[dim]; // initial state of point 0
    for(i = 0; i < dim; ++i)
		initVec[i] = i;

    nPts = ptsPerDim * dim;

    printf("Using up to %d threads\n", nthreads);

    State *state = State_new(nPts, dim);
    assert(state);
    for(i=0; i < nPts; ++i) {
        double x[dim];
        initState(i, initVec, dim, x, rng);
        State_setVector(state, i, dim, x);
    }

    // Check inputs
    if(dim <= 0) {
        fprintf(stderr,
                "%s:%d:Err dim=%d, should be > 0\n", __FILE__, __LINE__, dim);
        exit(1);
    }
    if(nPts <= 0) {
        fprintf(stderr,
                "%s:%d:Err nPts=%d, should be > 0\n", __FILE__, __LINE__,
                nPts);
        exit(1);
    }

    if((CR < 0) || (CR > 1.0)) {
        fprintf(stderr,
                "%s:%d:Err CR=%f, should be in [0,1]\n", __FILE__, __LINE__,
                CR);
        exit(1);
    }
    if(refresh <= 0) {
        fprintf(stderr,"%s:%d:Err refresh=%d, should be > 0\n",
                __FILE__, __LINE__, refresh);
        exit(1);
    }
    if(genmax <= 0) {
        fprintf(stderr,"%s:%d:Err genmax=%d, should be > 0\n",
                __FILE__, __LINE__, genmax);
        exit(1);
    }

    if(verbose) {
        printf("Strategy: %s\n", diffEvStrategyLbl(strategy));
        printf("nPts=%d F=%-4.2lg CR=%-4.2lg\n", nPts, F, CR);
    }

    // Number of parameters in quadratic model used to estimate
    // Hessian matrix: 1 intercept
    // dim linear terms
    // dim squared terms
    // (dim*(dim-1))/2 cross product terms
    unsigned nQuadPar = 1 + 2*dim + (dim*(dim-1))/2;

    // Number of points to use in fitting quadratic model
    unsigned nQuadPts = 10*nQuadPar;

    // parameters for Differential Evolution
    DiffEvPar   dep = {
        .dim = dim,
        .ptsPerDim = ptsPerDim,
        .refresh = refresh,
        .strategy = strategy,
        .nthreads = nthreads,
        .verbose = verbose,
        .seed = ((unsigned long) time(NULL))-1ul,
        .F = F,
        .CR = CR,
        .jobData = NULL,
        .JobData_dup = NULL,
        .JobData_free = NULL,
        .objfun = objFunc,
        .threadData = NULL,
        .ThreadState_new = NULL,
        .ThreadState_free = NULL,
		.state = state,
        .simSched = simSched,
        .ytol = ytol,
        .pb = PointBuff_new(dim, nQuadPts)
    };

    double      estimate[dim];
    double      cost, yspread;

    int         status = diffev(dim, estimate, &cost, &yspread, dep, rng);
    int ok=0;
    switch (status) {
    case 0:
        if(verbose) {
            printf("DiffEv converged. cost=%0.5lg yspread=%0.5lg\n",
                   cost, yspread);
            printf("Fitted parameters:");
            for(i = 0; i < dim; ++i)
                printf(" %lf", estimate[i]);
            putchar('\n');
        }
        ok = 1;
        break;
    default:
        ok = 0;
        break;
    }

    SimSched_free(simSched);
    gsl_rng_free(rng);

    unitTstResult("diffev", ok ? "OK" : "FAIL");

    return 0;
}
