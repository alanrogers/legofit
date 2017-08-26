/**
@file legofit.c
@page legofit
@brief Estimate parameters describing population sizes, the times
of separations and of episodes of gene flow, and levels of gene flow.

# `legofit`: estimate population history from site pattern data

    usage: legofit [options] input.lgo sitepat.txt
       where file input.lgo describes population history,
       and file sitepat.txt contains site pattern frequencies.
    Options may include:
       -T <x> or --tol <x>
          termination criterion
       -t <x> or --threads <x>
          number of threads (default is auto)
       -F <x> or --scaleFactor <x>
          set DE scale factor
       -x <x> or --crossover <x>
          set DE crossover probability
       -s <x> or --strategy <x>
          set DE strategy
       -S <g>@<r> or --stage <g>@<r>
          add stage with <g> generations and <r> simulation reps
       -p <x> or --ptsPerDim <x>
          number of DE points per free var
       -1 or --singletons
          Use singleton site patterns
       -v or --verbose
          verbose output
       -h or --help
          print this message

Two arguments are required:

 - an input file in @ref lgo ".lgo" format, which describes the history
   of population size, subdivision, and gene flow;
 - a file in the format produced by @ref tabpat "tabpat", which
   provides counts of site patterns in the data.

Legofit estimates the values of all the parameters that are declared
as "free" (rather than "fixed", "gaussian", or "constrained") in the
.lgo file. It does this by minimizing a "cost function", which
measures the difference between observed and expected
values. Currently, the cost function is the negative of composite
likelihood.  Other options are available via compile-time option (see
`typedefs.h`). I don't yet know which cost function is best.

Expected counts are estimated by computer simulation, and optimization
is done using the "differential evolution" (DE) algorithm.  The DE
algorithm maintains a swarm of points, each at a different set of
parameter values. The objective function is evaluated at these points
in a multithreaded job queue, so the program runs fasted on a machine
with lots of cores. You can set the number of threads using the `-t`
argument. By default, the program uses as many threads as there are
processors on the machine---usually the number of hypercores. The
optimal number of threads is usually somewhat smaller that this
default.

The DE algorithm can be tuned via command line arguments `-F`, `-x`,
`-s`, and `-p`. Details regarding these choices can be found in
"Differential evolution: a practical approach to global optimization",
by Price, Storn, and Lampinen. I don't yet know what is best, but I'm
currently using `-s 2` and `-s 4`, with default values of the other
options.

During the first few hundred generations of the DE algorithm, the
swarm of points adapts to the objective function. During this initial
phase, high accuracy is not required, so it can speed things up to use
low resolution in function evaluations. This is the purpose of the
`-S` argument, which adds a "stage" to the simulation schedule. For
example, `-S 1000@10000` adds a stage of 1000 DE generations, in each
of which 10000 simulation replicates are used to evaluate each
objective function. The `-S` argument can be given several times to
set up a simulation schedule with several stages. The algorithm is
allowed to converge only during the final stage. I am currently using
a 2-stage schedule: `-S 1000@10000 -S `1000@2000000`.

The `-1` option tells legofit to use singleton site patterns--patterns
in which the derived allele is present in only a single sample. This
is a bad idea with low-coverage sequence data. It also behaves poorly
with archaic DNA, because singletons are likely to be sequencing
artifacts. With high-quality modern DNA, it may one day prove
useful. To use it, you need to generate a data set that contains
singletons, by using the `-1` option of @\ref tabpat "tabpat".

2017-08-26: The convergence criterion has been changed. In the
previous code, differential evolution (DE) iterations stopped when
neither the best objective function value nor the spread of these
values had changed in a fixed number of iterations. That criterion was
used in our recent paper, "Early history of Neanderthals and
Denisovans", which was just published in PNAS.

I began to notice convergence problems with models larger than
those used in the PNAS paper. All bootstrap replicates would report
convergence, but some yielded wild parameter estimates. In these
outliers, the spread of objective function values was also very
large, indicating that the algorithm had not really converged. So I
implemented a new convergence criterion, based on the spread of
objective function values. The iterations terminate when this
spread falls to a pre-determined value.

This new convergence criterion works best with the KL
(Kullback-Leibler) cost function. Minimizing KL is the same as
maximizing log composite likelihood (lnL). But KL equals zero when
the fit is perfect. For this reason, it isn't necessary to scale
the tolerance value to lnL.

@copyright Copyright (c) 2016, 2017, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
 */

#include "branchtab.h"
#include "cost.h"
#include "diffev.h"
#include "gptree.h"
#include "lblndx.h"
#include "parstore.h"
#include "patprob.h"
#include "simsched.h"
#include <assert.h>
#include <float.h>
#include <getopt.h>
#include <gsl/gsl_sf_gamma.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

extern pthread_mutex_t seedLock;
extern unsigned long rngseed;
extern volatile sig_atomic_t sigstat;

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
#if COST==KL_COST || COST==LNL_COST
    fprintf(stderr,"usage: legofit [options] input.lgo sitepat.txt\n");
    fprintf(stderr,"   where file input.lgo describes population history,\n"
            "   and file sitepat.txt contains site pattern frequencies.\n");
#else
    fprintf(stderr,"usage: legofit [options] -u <mut_rate>"
            " -n <genome_size> input.lgo sitepat.txt\n");
    fprintf(stderr,"   where <mut_rate> is the mutation rate per nucleotide\n"
            "   site per generation, <genome_size> is the number of\n"
            "   nucleotides per haploid genome, file input.lgo describes\n"
            "   population history, and file sitepat.txt contains site\n"
            "   pattern frequencies.\n");
#endif
    fprintf(stderr,"Options may include:\n");
    tellopt("-T <x> or --tol <x>", "termination criterion");
    tellopt("-t <x> or --threads <x>", "number of threads (default is auto)");
    tellopt("-F <x> or --scaleFactor <x>", "set DE scale factor");
    tellopt("-x <x> or --crossover <x>", "set DE crossover probability");
    tellopt("-s <x> or --strategy <x>", "set DE strategy");
    tellopt("-S <g>@<r> or --stage <g>@<r>",
            "add stage with <g> generations and <r> simulation reps");
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

    // Install handler for keyboard interrupts.
    signal(SIGINT, sighandle);

    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
        {"threads", required_argument, 0, 't'},
		{"crossover", required_argument, 0, 'x'},
		{"scaleFactor", required_argument, 0, 'F'},
        {"strategy", required_argument, 0, 's'},
        {"stage", required_argument, 0, 'S'},
        {"tol", required_argument, 0, 'T'},
        {"ptsPerDim", required_argument, 0, 'p'},
#if COST!=KL_COST && COST!=LNL_COST
        {"mutRate", required_argument, 0, 'u'},
        {"genomeSize", required_argument, 0, 'n'},
#endif
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
    double      lo_twoN = 1.0, hi_twoN = 1e7;  // twoN bounds
    double      lo_t = 0.0, hi_t = 1e7;        // t bounds
    int         nThreads = 0;     // total number of threads
    int         doSing=0;  // nonzero means use singleton site patterns
    int         status, optndx;
    long        simreps = 1000000;
    char        lgofname[200] = { '\0' };
    char        patfname[200] = { '\0' };

	// DiffEv parameters
	double      F = 0.9;
	double      CR = 0.8;
#if COST!=KL_COST && COST!=LNL_COST
    double      u = 0.0;       // mutation rate per site per generation
    long        nnuc = 0;      // number of nucleotides per haploid genome
#endif
    double      ytol = 1e-4;   // stop when yspread <= ytol
	int         strategy = 1;
	int         ptsPerDim = 10;
    int         verbose = 0;
    SimSched    *simSched = SimSched_new();

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
#if COST==KL_COST || COST==LNL_COST
        i = getopt_long(argc, argv, "T:t:F:p:s:S:a:vx:1h",
                        myopts, &optndx);
#else
        i = getopt_long(argc, argv, "T:t:F:p:s:S:a:vx:u:n:1h",
                        myopts, &optndx);
#endif
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            usage();
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
		case 'S':
            {
                // Add a stage to simSched.
                char b[20], *g, *r;
                status = snprintf(b, sizeof b, "%s", optarg);
                if(status >= sizeof b) {
                    fprintf(stderr,"%s:%d: buffer overflow reading arg %s\n",
                            __FILE__,__LINE__,optarg);
                    exit(EXIT_FAILURE);
                }
                g = r = b;
                (void) strsep(&r, "@");
                if(r==NULL
                   || strlen(r) == 0
                   || strlen(g) == 0)
                    usage();
                long stageGen = strtol(g, NULL, 10);
                long stageRep  = strtol(r, NULL, 10);
                simreps = stageRep;
                SimSched_append(simSched, stageGen, stageRep);
            }
            break;
        case 'v':
            verbose = 1;
            break;
        case 'T':
            ytol = strtod(optarg, 0);
            break;
		case 'x':
			CR = strtod(optarg, 0);
			break;
#if COST!=KL_COST && COST!=LNL_COST
        case 'u':
            u = strtod(optarg, 0);
            break;
        case 'n':
            nnuc = strtol(optarg, NULL, 10);
            break;
#endif
        case '1':
            doSing=1;
            break;
        case 'h':
            usage();
            break;
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

#if COST!=KL_COST && COST!=LNL_COST
    if(u==0.0) {
        fprintf(stderr,"Use -u to set mutation rate per generation.\n");
        usage();
    }

    if(nnuc==0) {
        fprintf(stderr,"Use -n to set # of nucleotides per haploid genome.\n");
        usage();
    }
#endif

    snprintf(lgofname, sizeof(lgofname), "%s", argv[optind]);
    assert(lgofname[0] != '\0');
    snprintf(patfname, sizeof(patfname), "%s", argv[optind+1]);
    assert(patfname[0] != '\0');

    // Default simulation schedule.
    // Stage 1: 200 DE generations of 1000 simulation replicates
    // Stage 2: 100 generations of 10000 replicates
    if(0 == SimSched_nStages(simSched)) {
        SimSched_append(simSched, 200, 1000);
        SimSched_append(simSched, 100, 10000);
        SimSched_append(simSched, 1000, simreps);
    }

    SimSched_print(simSched, stdout);

    Bounds bnd = {
            .lo_twoN = lo_twoN,
            .hi_twoN = hi_twoN,
            .lo_t = lo_t,
            .hi_t = hi_t
    };

    GPTree *gptree = GPTree_new(lgofname, bnd);
	LblNdx lblndx  = GPTree_getLblNdx(gptree);

    int dim = GPTree_nFree(gptree); // number of free parameters
    if(dim == 0) {
        fprintf(stderr,"Error@%s:%d: no free parameters\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    if(nThreads == 0)
        nThreads = ceil(0.75*getNumCores());
    if(nThreads > dim*ptsPerDim)
        nThreads = dim*ptsPerDim;

    printf("# DE strategy        : %d\n", strategy);
    printf("#    F               : %lf\n", F);
    printf("#    CR              : %lf\n", CR);
    printf("#    tolerance       : %lf\n", ytol);
    printf("# nthreads           : %d\n", nThreads);
    printf("# lgo input file     : %s\n", lgofname);
    printf("# site pat input file: %s\n", patfname);
    printf("# pts/dimension      : %d\n", ptsPerDim);
#if COST!=KL_COST && COST!=LNL_COST
    printf("# mut_rate/generation: %lg\n", u);
    printf("# nucleotides/genome : %ld\n", nnuc);
#endif
    printf("# %s singleton site patterns.\n",
           (doSing ? "Including" : "Excluding"));
#if COST==KL_COST
    printf("# cost function      : %s\n", "KL");
#elif COST==LNL_COST
    printf("# cost function      : %s\n", "negLnL");
#elif COST==CHISQR_COST
    printf("# cost function      : %s\n", "ChiSqr");
#elif COST==SMPLCHISQR_COST
    printf("# cost function      : %s\n", "SmplChiSqr");
#elif COST==POISSON_COST
    printf("# cost function      : %s\n", "Poisson");
#else
# error "Unknown cost function"
#endif

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
        .doSing = doSing,
#if COST!=KL_COST && COST!=LNL_COST
        .u = u,
        .nnuc = nnuc,
#endif
        .simSched = simSched
    };

    // parameters for Differential Evolution
    DiffEvPar   dep = {
        .dim = dim,
        .ptsPerDim = ptsPerDim,
        .refresh = 2,  // how often to print a line of output
        .strategy = strategy,
        .nthreads = nThreads,
        .verbose = verbose,
        .seed = ((unsigned long) time(NULL))-1ul,
        .F = F,
        .CR = CR,
		.jobData = &costPar,
        .JobData_dup = CostPar_dup,
        .JobData_free = CostPar_free,
        .objfun = costFun,
		.threadData = NULL,
		.ThreadState_new = ThreadState_new,
		.ThreadState_free = ThreadState_free,
        .initData = gptree,
        .initialize = initStateVec,
        .simSched = simSched,
        .ytol = ytol
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

    status = diffev(dim, estimate, &cost, &yspread, dep, rng);

    printf("DiffEv %s. cost=%0.5lg spread=%0.5lg\n",
           status==0 ? "converged" : "FAILED", cost, yspread);
#if COST==LNL_COST
    printf("  relspread=%e", yspread/cost);
#endif
    putchar('\n');

    // Get mean site pattern branch lengths
    GPTree_setParams(gptree, dim, estimate);
    BranchTab *bt = patprob(gptree, simreps, doSing, rng);
    BranchTab_divideBy(bt, (double) simreps);
    //    BranchTab_print(bt, stdout);

    printf("Fitted parameter values\n");
#if 1
	GPTree_printParStoreFree(gptree, stdout);
#else
	GPTree_printParStore(gptree, stdout);
#endif

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
    BranchTab_free(obs);
    gsl_rng_free(rng);
    GPTree_sanityCheck(gptree, __FILE__, __LINE__);
    GPTree_free(gptree);
    SimSched_free(simSched);
    fprintf(stderr,"legofit is finished\n");

    return 0;
}
