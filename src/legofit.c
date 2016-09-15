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
#include <float.h>
#include <pthread.h>
#include <gsl/gsl_sf_gamma.h>

extern pthread_mutex_t seedLock;
extern unsigned long rngseed;

void        usage(void);
void        initStateVec(int ndx, void *void_p, int n, double x[n],
                         gsl_rng *rng);
void       *ThreadState_new(void *notused);
void        ThreadState_free(void *rng);

void *ThreadState_new(void *notused) {
	// Lock seed, initialize random number generator, increment seed,
	// and unlock.
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);

	pthread_mutex_lock(&seedLock);
    gsl_rng_set(rng, rngseed);
	rngseed = (rngseed == ULONG_MAX ? 0 : rngseed+1);
	pthread_mutex_unlock(&seedLock);

    return rng;
}

void ThreadState_free(void *rng) {
    gsl_rng_free( (gsl_rng *) rng );
}

void usage(void) {
    fprintf(stderr,"usage: legofit [options] -u <mut_rate>"
            " -n <genome_size> input.lgo sitepat.txt\n");
    fprintf(stderr,"   where <mut_rate> is the mutation rate per nucleotide\n"
            "   site per generation, <genome_size> is the number of\n"
            "   nucleotides per haploid genome, file input.lgo describes\n"
            "   population history, and file sitepat.txt contains site\n"
            "   pattern frequencies.\n");
    fprintf(stderr,"Options may include:\n");
    tellopt("-i <x> or --deItr <x>", "number of DE iterations");
    tellopt("-r <x> or --simreps <x>", "number of reps in each function eval");
    tellopt("-a <x> or --DEtolerance <x>", "termination criterion");
    tellopt("-t <x> or --threads <x>", "number of threads (default is auto)");
    tellopt("-F <x> or --scaleFactor <x>", "set DE scale factor");
    tellopt("-x <x> or --crossover <x>", "set DE crossover probability");
    tellopt("-s <x> or --strategy <x>", "set DE strategy");
    tellopt("-p <x> or --ptsPerDim <x>", "number of DE points per free var");
	tellopt("-1 or --singletons", "Use singleton site patterns");
    tellopt("-v or --verbose", "verbose output");
    tellopt("-h or --help", "print this message");
    exit(1);
}

/// Initialize vector x. If ndx==0, simply copy the parameter vector
/// from the GPTree object. Otherwise, randomize the GPTree first.
/// This ensures that differential evolution starts with a set of
/// points, one of which is the same as the values in the input
/// file. This allows you to improve on existing estimates without
/// starting from scratch each time.
void initStateVec(int ndx, void *void_p, int n, double x[n], gsl_rng *rng){
    GPTree *gpt = (GPTree *) void_p;
    if(ndx == 0)
        GPTree_getParams(gpt, n, x);
    else {
        GPTree *g2 = GPTree_dup(gpt);
        GPTree_randomize(g2, rng);
        GPTree_getParams(g2, n, x);
        GPTree_free(g2);
    }
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
        {"DEtolerance", required_argument, 0, 'a'},
        {"ptsPerDim", required_argument, 0, 'p'},
        {"mutRate", required_argument, 0, 'u'},
        {"genomeSize", required_argument, 0, 'n'},
        {"singletons", no_argument, 0, '1'},
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, 0, 'v'},
        {NULL, 0, NULL, 0}
    };

    printf("########################################\n"
           "# legofit: estimate population history #\n"
           "########################################\n");
    putchar('\n');

    int         i, j;
    time_t      currtime = time(NULL);
	unsigned long pid = (unsigned long) getpid();
    double      lo_twoN = 0.0, hi_twoN = 1e6;  // twoN bounds
    double      lo_t = 0.0, hi_t = 1e6;        // t bounds
    int         nThreads = 0;     // total number of threads
    int         doSing=0;  // nonzero means use singleton site patterns
    int         optndx;
    long        simreps = 100;
    char        lgofname[200] = { '\0' };
    char        patfname[200] = { '\0' };

	// DiffEv parameters
	double      F = 0.9;
	double      CR = 0.8;
	double      DEtol = 0.05; // termination criterion
    double      u = 0.0;      // mutation rate per site per generation
    long        nnuc = 0;     // number of nucleotides per haploid genome
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
        i = getopt_long(argc, argv, "i:t:F:p:r:s:a:vx:u:n:1h", myopts, &optndx);
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
            DEtol = strtod(optarg, 0);
            break;
		case 'x':
			CR = strtod(optarg, 0);
			break;
        case 'u':
            u = strtod(optarg, 0);
            break;
        case 'n':
            nnuc = strtol(optarg, NULL, 10);
            break;
        case '1':
            doSing=1;
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

    if(u==0.0) {
        fprintf(stderr,"Use -u to set mutation rate per generation.\n");
        usage();
    }

    if(nnuc==0) {
        fprintf(stderr,"Use -n to set # of nucleotides per haploid genome.\n");
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
    printf("#    deTol           : %lf\n", DEtol);
    printf("#    F               : %lf\n", F);
    printf("#    CR              : %lf\n", CR);
    printf("# simreps            : %lu\n", simreps);
    printf("# nthreads           : %d\n", nThreads);
    printf("# lgo input file     : %s\n", lgofname);
    printf("# site pat input file: %s\n", patfname);
    printf("# pts/dimension      : %d\n", ptsPerDim);
    printf("# mut_rate/generation: %lg\n", u);
    printf("# nucleotides/genome : %ld\n", nnuc);
    printf("# %s singleton site patterns.\n",
           (doSing ? "Including" : "Excluding"));
#if COST==KL_COST
    printf("# cost function      : %s\n", "KL");
#elif COST==CHISQR_COST
    printf("# cost function      : %s\n", "ChiSqr");
#elif COST==SMPLCHISQR_COST
    printf("# cost function      : %s\n", "SmplChiSqr");
#elif COST==POISSON_COST
    printf("# cost function      : %s\n", "Poisson");
#else
# error "Unknown cost function"
#endif

    Bounds bnd = {
            .lo_twoN = lo_twoN,
            .hi_twoN = hi_twoN,
            .lo_t = lo_t,
            .hi_t = hi_t
    };
    GPTree *gptree = GPTree_new(lgofname, bnd);
	LblNdx lblndx  = GPTree_getLblNdx(gptree);

    // Observed site pattern frequencies
    BranchTab *obs = BranchTab_parse(patfname, &lblndx);
    if(doSing) {
        if(!BranchTab_hasSingletons(obs)) {
            fprintf(stderr,"%s:%d: Command line includes singletons "
                    "(-1 or --singletons)\n"
                    "    but none are present in \"%s\".\n",
                    __FILE__,__LINE__, patfname);
            exit(EXIT_FAILURE);
        }
    }else{
        if(BranchTab_hasSingletons(obs)) {
            fprintf(stderr,"%s:%d: Command line excludes singletons "
                    "(neither -1 nor --singletons)\n"
                    "    but singletons are present in \"%s\".\n",
                    __FILE__,__LINE__, patfname);
            exit(EXIT_FAILURE);
        }
    }
#if COST==KL_COST
    BranchTab_normalize(obs);
#endif

    // parameters for cost function
    CostPar costPar = {
        .obs = obs,
        .gptree = gptree,
        .nThreads = nThreads,
        .nreps = simreps,
        .doSing = doSing,
        .u = u,
        .nnuc = nnuc
    };

    // parameters for Differential Evolution
    int dim = GPTree_nFree(gptree); // number of free parameters
    if(dim == 0) {
        fprintf(stderr,"Error@%s:%d: no free parameters\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    DiffEvPar   dep = {
        .dim = dim,
        .ptsPerDim = ptsPerDim,
        .genmax = deItr,
        .refresh = 2,  // how often to print a line of output
        .strategy = strategy,
        .nthreads = nThreads,
        .verbose = verbose,
        .seed = ((unsigned long) time(NULL))-1ul,
        .F = F,
        .CR = CR,
        .DEtol = DEtol,
		.jobData = &costPar,
        .JobData_dup = CostPar_dup,
        .JobData_free = CostPar_free,
        .objfun = costFun,
		.threadData = NULL,
		.ThreadState_new = ThreadState_new,
		.ThreadState_free = ThreadState_free,
        .initData = gptree,
        .initialize = initStateVec
    };

    double      estimate[dim];
    double      cost, yspread;

    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, rngseed);
	rngseed = (rngseed == ULONG_MAX ? 0 : rngseed+1);

    printf("Initial parameter values\n");
	GPTree_printParStore(gptree, stdout);

    // Flush just before diffev so output file will be as complete as
    // possible while diffev is running.
    fflush(stdout);

    int         status = diffev(dim, estimate, &cost, &yspread, dep, rng);
    switch (status) {
    case 0:
        printf("DiffEv Converged. cost=%0.5lg; spread=%0.5lg < %lg = deTol\n",
               cost, yspread, DEtol);
        break;
    default:
        printf("DiffEv FAILED. cost=%0.5lg; spread=%0.5lg > %lg = deTol\n",
               cost, yspread, DEtol);
        break;
    }

    // Get mean site pattern branch lengths
    GPTree_setParams(gptree, dim, estimate);
    BranchTab *bt = patprob(gptree, nThreads, simreps, doSing, rng);
    BranchTab_divideBy(bt, simreps);
    
    printf("Fitted parameter values\n");
	GPTree_printParStoreFree(gptree, stdout);

    // Put site patterns and branch lengths into arrays.
    unsigned npat = BranchTab_size(bt);
    tipId_t pat[npat];
    double brlen[npat];
    double sqr[npat];
    BranchTab_toArrays(bt, npat, pat, brlen, sqr);

    // Determine order for printing lines of output
    unsigned ord[npat];
    orderpat(npat, ord, pat);

    printf("#%14s %10s\n", "SitePat", "BranchLen");
    char        buff[100];
    for(j = 0; j < npat; ++j) {
        char        buff2[100];
        snprintf(buff2, sizeof(buff2), "%s",
                 patLbl(sizeof(buff), buff, pat[ord[j]], &lblndx));
        printf("%15s %10.7lf\n", buff2, brlen[ord[j]]);
    }

    BranchTab_free(bt);
    gsl_rng_free(rng);
    GPTree_free(gptree);
    fprintf(stderr,"legofit is finished\n");

    return 0;
}
