#include "diffev.h"
#include "misc.h"
#include <getopt.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <memory.h>
#include <gsl/gsl_rng.h>

void        usage(void);
double      objFunc(int dim, double x[dim], void *jdat, void *tdat);

#define RUGGED

/// Local minima wherever all entries of v are integers. Global
/// minimum at v=0.
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

#undef RUGGED

/// Print usage message and exit.
void usage(void) {
    fprintf(stderr, "usage: diffev [options]\n");
    fprintf(stderr, "   where options may include:\n");

    /* misc */
    tellopt("-s <method> or --strategy <method>", "strategy");
    tellopt("-g <x> or --genmax <x>", "max generations");
    tellopt("-r <x> or --refresh <x>", "refresh interval");
    tellopt("-n or --nParam", " number of parameters");
    tellopt("-p <x> or --ptsPerDim <x>", "points per dimension");
    tellopt("-F <x> or --F <x>", "DE weight factor");
    tellopt("-c <x> or --crossOver <x>", "crossover probability");
    tellopt("-t <x> or --threads <x>", "number of threads (default is auto)");
    tellopt("-v or --verbose", "more output");
    tellopt("-h or --help", "print this message");
    exit(1);
}

int main(int argc, char *argv[]) {
    int         optndx;
    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
        {"strategy", required_argument, 0, 's'},
        {"threads", required_argument, 0, 't'},
        {"genmax", required_argument, 0, 'g'},
        {"refresh", required_argument, 0, 'r'},
        {"nParam", required_argument, 0, 'n'},
        {"ptsPerDim", required_argument, 0, 'p'},
        {"F", required_argument, 0, 'F'},
        {"crossOver", required_argument, 0, 'c'},
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
    double      deTol = 1e-4;   // tolerance

    int         i;
    int         nthreads = 0;
    time_t      currtime = time(NULL);
    unsigned    baseSeed = currtime % UINT_MAX;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    CHECKMEM(rng);
    gsl_rng_set(rng, baseSeed);

    // command line arguments
    for(;;) {
        i = getopt_long(argc, argv, "s:t:g:r:n:p:F:c:hv", myopts, &optndx);
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
        case 't':
            nthreads = strtol(optarg, NULL, 10);
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

    // There should be no non-flag arguments
    if(argc > optind)
        usage();

    // lower and upper bounds on each parameter
    double      loBound[dim], hiBound[dim];
    for(i = 0; i < dim; ++i) {
        loBound[i] = -10.0;
        hiBound[i] = 10.0;
    }

    nPts = ptsPerDim * dim;

    if(nthreads == 0)
        nthreads = getNumCores();

    printf("Using up to %d threads\n", nthreads);

    // Check inputs
    if(dim <= 0)
        eprintf("%s:%d:Err dim=%d, should be > 0\n", __FILE__, __LINE__, dim);
    if(nPts <= 0)
        eprintf("%s:%d:Err nPts=%d, should be > 0\n", __FILE__, __LINE__,
                nPts);

    if((CR < 0) || (CR > 1.0))
        eprintf("%s:%d:Err CR=%f, should be in [0,1]\n", __FILE__, __LINE__,
                CR);
    if(refresh <= 0)
        eprintf("%s:%d:Err refresh=%d, should be > 0\n",
                __FILE__, __LINE__, refresh);
    if(genmax <= 0)
        eprintf("%s:%d:Err genmax=%d, should be > 0\n",
                __FILE__, __LINE__, genmax);
    for(i = 0; i < dim; ++i) {
        if(loBound[i] > hiBound[i])
            eprintf("%s:%d:Err loBound[%d]=%lf > hiBound[%d]=%lf\n",
                    __FILE__, __LINE__, i, loBound[i], i, hiBound[i]);
    }

    printf("Strategy: %s\n", diffEvStrategyLbl(strategy));
    printf("nPts=%d F=%-4.2lg CR=%-4.2lg\n", nPts, F, CR);

    // parameters for Differential Evolution
    DiffEvPar   dep = {
        .dim = dim,
        .ptsPerDim = ptsPerDim,
        .genmax = genmax,
        .refresh = refresh,
        .strategy = strategy,
        .nthreads = nthreads,
        .seed = (unsigned long) time(NULL),
        .verbose = verbose,
        .F = F,
        .CR = CR,
        .deTol = deTol,
        .loBound = loBound,
        .hiBound = hiBound,
        .jobData = NULL,
        .objfun = objFunc,
        .threadData = NULL,
        .ThreadState_new = NULL,
        .ThreadState_free = NULL
    };

    double      estimate[dim];
    double      cost, yspread;

    int         status = diffev(dim, estimate, &cost, &yspread, dep, rng);
    switch (status) {
    case 0:
        printf("DiffEv converged. cost=%0.5lg costSpread=%0.5lg\n",
               cost, yspread);
        printf("Fitted parameters:");
        for(i = 0; i < dim; ++i)
            printf(" %lf", estimate[i]);
        putchar('\n');
        break;
    default:
        printf("DiffEv FAILED\n");
        break;
    }

    gsl_rng_free(rng);

    return 0;
}
