/**
 * @file legofit.c
 * @brief Estimate parameters describing population sizes, the times of separations
 * and of episodes of gene flow, and levels of gene flow.
 *
 * @copyright Copyright (c) 2016, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "branchtab.h"
#include "cost.h"
#include "diffev.h"
#include "gptree.h"
#include "lblndx.h"
#include "parstore.h"
#include "patprob.h"
#include <assert.h>
#include <getopt.h>
#include <limits.h>
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

extern unsigned long rngseed;

void        usage(void);

void usage(void) {
    fprintf(stderr,"usage: legofit [options] input.lgo sitepat.txt\n");
    fprintf(stderr,"   where file input.lgo describes population history,\n");
    fprintf(stderr,"   file sitepat.txt contains site pattern frequencies,\n");
    fprintf(stderr,"   and options may include:\n");
    tellopt("-i <x> or --deItr <x>", "number of DE iterations");
    tellopt("-r <x> or --simreps <x>", "number of reps in each function eval");
    tellopt("-a <x> or --deTol <x>", "DE tolerance: smaller means less accurate");
    tellopt("-t <x> or --threads <x>", "number of threads (default is auto)");
    tellopt("-F <x> or --scaleFactor <x>", "set DE scale factor");
    tellopt("-x <x> or --crossover <x>", "set DE crossover probability");
    tellopt("-v or --verbose", "verbose output");
    tellopt("-h or --help", "print this message");
    exit(1);
}

int main(int argc, char **argv) {

    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
        {"deItr", required_argument, 0, 'i'},
        {"threads", required_argument, 0, 't'},
		{"crossover", required_argument, 0, 'x'},
		{"scaleFactor", required_argument, 0, 'F'},
		{"simreps", required_argument, 0, 'r'},
        {"strategy", required_argument, 0, 's'},
        {"deTol", required_argument, 0, 'a'},
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, 0, 'v'},
        {NULL, 0, NULL, 0}
    };

    printf("########################################\n"
           "# legofit: estimate population history #\n"
           "########################################\n");
    putchar('\n');

    int         i;
    time_t      currtime = time(NULL);
	unsigned long pid = (unsigned long) getpid();
    double      lo_twoN = 0.0, hi_twoN = 1e6;  // twoN bounds
    double      lo_t = 0.0, hi_t = 1e6;        // t bounds
    int         nThreads = 0;     // total number of threads
    int         optndx;
    long        simreps = 100;
    char        lgofname[200] = { '\0' };
    char        patfname[200] = { '\0' };

	// DiffEv parameters
	double      F = 0.9;
	double      CR = 0.8;
	double      deTol = 1e-3;
    int         deItr = 1000; // number of diffev iterations
	int         strategy = 1;
	int         ptsPerDim = 10;
    int         verbose = 0;

#if defined(__DATE__) && defined(__TIME__)
    printf("# Program was compiled: %s %s\n", __DATE__, __TIME__);
#endif
    printf("# Program was run: %s\n", ctime(&currtime));

    printf("# cmd:");
    for(i = 0; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');
    fflush(stdout);

	rngseed = currtime^pid;

    // command line arguments
    for(;;) {
        i = getopt_long(argc, argv, "i:t:F:p:r:s:a:vx:h", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            usage();
            break;
        case 'r':
            simreps = strtol(optarg, 0, 10);
            break;
        case 'i':
            deItr = strtol(optarg, 0, 10);
            break;
        case 't':
            nThreads = strtol(optarg, NULL, 10);
            break;
		case 'F':
			F = strtod(optarg, 0);
			break;
        case 'p':
            ptsPerDim = strtol(optarg, NULL, 10);
            break;
		case 's':
			strategy = strtol(optarg, NULL, 10);
            break;
        case 'v':
            verbose = 1;
            break;
        case 'a':
            deTol = strtod(optarg, 0);
            break;
		case 'x':
			CR = strtod(optarg, 0);
			break;
        case 'h':
        default:
            fprintf(stderr,"Can't parse option %c\n", i);
            usage();
        }
    }

    // remaining options gives file names
    if(argc - optind != 2) {
        fprintf(stderr, "Command line must specify 2 input files.\n");
        usage();
    }
        
    snprintf(lgofname, sizeof(lgofname), "%s", argv[optind]);
    assert(lgofname[0] != '\0');
    snprintf(patfname, sizeof(patfname), "%s", argv[optind+1]);
    assert(patfname[0] != '\0');

    if(nThreads == 0)
        nThreads = getNumCores();

    if(nThreads > simreps)
        nThreads = simreps;

    printf("# DE strategy        : %d\n", strategy);
    printf("#    deItr           : %d\n", deItr);
    printf("#    deTol           : %lf\n", deTol);
    printf("#    F               : %lf\n", F);
    printf("#    CR              : %lf\n", CR);
    printf("# simreps            : %lu\n", simreps);
    printf("# nthreads           : %d\n", nThreads);
    printf("# lgo input file     : %s\n", lgofname);
    printf("# site pat input file: %s\n", patfname);

    Bounds bnd = {
            .lo_twoN = lo_twoN,
            .hi_twoN = hi_twoN,
            .lo_t = lo_t,
            .hi_t = hi_t
    };
    GPTree *gptree = GPTree_new(lgofname, bnd);
	LblNdx lblndx  = GPTree_getLblNdx(gptree);

	GPTree_printParStore(gptree, stdout);

    // Observed site pattern frequencies
    BranchTab *obs = BranchTab_parse(patfname, &lblndx);
    BranchTab_normalize(obs);

    // parameters for cost function
    CostPar costPar = {
        .obs = obs,
        .gptree = gptree,
        .nThreads = nThreads,
        .nreps = simreps
    };

    // parameters for Differential Evolution
    int dim = GPTree_nFree(gptree); // number of free parameters
    DiffEvPar   dep = {
        .dim = dim,
        .ptsPerDim = ptsPerDim,
        .genmax = deItr,
        .refresh = 1,  // how often to print a line of output
        .strategy = strategy,
        .nthreads = 1,
        .verbose = verbose,
        .seed = ((unsigned long) time(NULL))-1ul,
        .F = F,
        .CR = CR,
        .deTol = deTol,
        .loBound = GPTree_loBounds(gptree),
        .hiBound = GPTree_upBounds(gptree),
		.jobData = &costPar,
        .objfun = costFun,
		.threadData = NULL,
		.ThreadState_new = NULL,
		.ThreadState_free = NULL,
        .randomizeData = gptree,
        .randomize = GPTree_randomize
    };

    double      estimate[dim];
    double      cost, yspread;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);

    gsl_rng_set(rng, rngseed);
	rngseed = (rngseed == ULONG_MAX ? 0 : rngseed+1);

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
    GPTree_free(gptree);

    return 0;
}
