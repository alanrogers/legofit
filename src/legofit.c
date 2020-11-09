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
       -e or --estimate
          estimate site pattern probabilities by simulation
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
       --stateIn <filename>
          read initial state of optimizer from file. Option may be repeated.
       --stateOut <filename>
          write final state of optimizer to file
       -1 or --singletons
          Use singleton site patterns
       -v or --verbose
          verbose output
       --version
          Print version and exit
       -h or --help
          print this message

Two arguments are required:

 - an input file in @ref lgo ".lgo" format, which describes the history
   of population size, subdivision, and gene flow;
 - a file in the format produced by @ref tabpat "tabpat", which
   provides counts of site patterns in the data.

Legofit estimates the values of all the parameters that are declared
as "free" in the .lgo file. It does this by minimizing a "cost
function", which measures the difference between observed and expected
values. Currently, the cost function is the negative of composite
likelihood.  Other options are available via compile-time option (see
`typedefs.h`). I don't yet know which cost function is best.

Expected counts are estimated by computer simulation, and optimization
is done using the "differential evolution" (DE) algorithm.  The DE
algorithm maintains a swarm of points, each at a different set of
parameter values. The objective function is evaluated at these points
in a multithreaded job queue, so the program runs fastest on a machine
with lots of cores. You can set the number of threads using the `-t`
argument. By default, the program uses 3/4 as many threads as there
are processors on the machine---usually the number of hypercores.

The DE algorithm can be tuned via command line arguments `-F`, `-x`,
`-s`, and `-p`. Details regarding these choices can be found in
"Differential evolution: a practical approach to global optimization",
by Price, Storn, and Lampinen. We've had good results with `-s 2`, `-F
0.3`, and `-x 0.8`, so these became the defaults as of version 1.15.

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

By default, the initial swarm of points consists of one point
representing the parameter values in the .lgo file, plus other points
scattered randomly throughout the feasible region of parameter
space. The total number of points defaults to 10 times the number of free
parameters. To change this number, see the --ptsPerDim option.

The initial swarm of points can also be specified using the
`--stateIn` option. This reads a file specifying the initial state of
the swarm of points maintained by DE. The format of this file is as
described below for the `--stateOut` option. The number of free
parameters in each `--stateIn` file should be as specified in the .lgo
file. The number of points in the files may differ.

The `--stateIn` option may be given more than once, each time with a
different input file. When more than one file is given, Legofit
constructs the initial swarm of points by combining points from all
input files.

The option `--stateOut` is used to define an output file for the final
state of the optimizer. This output file begins with a row giving the
number of points and the number of free parameters. After that, there
is a row for each point in the swarm of points maintained by
diffev.c. In each row, the first entry is the value of the cost
function at that point. The remaining entries give the free parameter
values in the same order in which they are printed by legofit.

The format of the state file changed in early July, 2018. Before that
date, the state file did not include names of parameters. Parameter
values in the state file had to be arranged in the same order as the
free parameters in the .lgo file. If the .lgo and state files referred
to different parameters or to the same parameters in a different
order, parameters would get the wrong values. No error was detected
unless this misassignment resulted in a tree that was not
feasible--for example, one in which a segment of the population tree
was older than its parent in the tree.

The new state file format includes the names of parameters, and the
input routine compares these against the names of free parameters in
the .lgo file. The two lists must have the same parameters in the same
order. Otherwise, legofit aborts with an error message. Old- and
new-format state files can both be input using --stateIn arguments and
can be intermingled in a single legofit run.

The `-1` option tells legofit to use singleton site patterns--patterns
in which the derived allele is present in only a single sample. This
is a bad idea with low-coverage sequence data, because singletons are
likely to be sequencing errors. We use this option routinely now with
high-coverage data. To use it, you need to generate a data set that
contains singletons, by using the `-1` option of @\ref tabpat
"tabpat", @ref sitepat "sitepat", @ref scrmpat "scrmpat", or @ref
simpat "simpat".

2017-08-26: The convergence criterion has been changed. In the
previous code, differential evolution (DE) iterations stopped when
neither the best objective function value nor the spread of these
values had changed in a fixed number of iterations. That criterion was
used in our PNAS paper, "Early history of Neanderthals and
Denisovans".

I began to notice convergence problems with models larger than those
used in the August 2017 PNAS paper. All bootstrap replicates would
report convergence, but some yielded wild parameter estimates. In
these outliers, the spread of objective function values was also very
large, indicating that the algorithm had not really converged. So I
implemented a new convergence criterion, based on the spread of
objective function values. The iterations terminate when this spread
falls to a pre-determined value.

This new convergence criterion works best with the KL
(Kullback-Leibler) cost function. Minimizing KL is the same as
maximizing log composite likelihood (lnL). But KL equals zero when
the fit is perfect. For this reason, it isn't necessary to scale
the tolerance value to lnL.

In each DE generation, the algorithm evaluates the cost function of
each point in the swarm. The cost function measures the difference
between observed and expected site pattern frequencies. The difference
between the maximum and minimum cost is called the "spread". The
algorithm stops when "spread" falls to the tolerance value,
"ytol". This tolerance can be set on the command line with "-T" or
"--tol" arguments. Convergence is checked only in the last stage, as
set by "-S" or "--stage".

If the algorithm fails to converge, there are several options. First,
you can use "-S" or "--stage" to increase the maximum number of DE
generations. For example, if the last stage is "-S 2000@2000000", you
are allowing 2000 DE generations, in each of which each cost function
value is estimated with 2000000 iterations. To double the allowed
number of DE generations, change this to "-S 4000@2000000".

Second, you can relax the tolerance. By default, this is 1e-4. It is
reported in the legofit output. To double this value, use "-T 2e-4" or
"--tol 2e-4".

Legofit handles three types of signal: SIGINT, SIGTERM, and
SIGUSR1. If legofit is running in the foreground, the first of these
signals can be generated by typing Ctrl-C. Otherwise (under linux or
osx), use `killall -SIGINT legofit`, `killall -SIGTERM legofit`, or
`killall -SIGUSR1 legofit`. In response to SIGINT or SIGTERM, Legofit
will wait until the end of the current generation of the diffev
algorithm and then exit gracefully, printing all the usual output.  In
response to SIGUSR1, Legofit will wait until the end of the current DE
generation and then write to stderr a summary of the state of the
optimizer.

@copyright Copyright (c) 2016, 2017, 2018, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
**/

#include "branchtab.h"
#include "comb.h"
#include "cost.h"
#include "diffev.h"
#include "lblndx.h"
#include "matcoal.h"
#include "network.h"
#include "parseopf.h"
#include "parstore.h"
#include "patprob.h"
#include "pointbuff.h"
#include "simsched.h"
#include "state.h"
#include <assert.h>
#include <float.h>
#include <getopt.h>
#include <gsl/gsl_sf_gamma.h>
#include <libgen.h>
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

void usage(void);
void *ThreadState_new(void *notused);
void ThreadState_free(void *rng);

void *ThreadState_new(void *notused) {
    // Lock seed, initialize random number generator, increment seed,
    // and unlock.
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);

    pthread_mutex_lock(&seedLock);
    gsl_rng_set(rng, rngseed);
    rngseed += 1; // wraps to 0 at ULONG_MAX
    pthread_mutex_unlock(&seedLock);

    return rng;
}

void ThreadState_free(void *rng) {
    gsl_rng_free((gsl_rng *) rng);
}

void usage(void) {
    fprintf(stderr, "usage: legofit [options] input.lgo sitepat.txt\n");
    fprintf(stderr, "   where file input.lgo describes population history,\n"
            "   and file sitepat.txt contains site pattern frequencies.\n");
    fprintf(stderr, "Options may include:\n");
    tellopt("-1 or --singletons", "Use singleton site patterns");
    tellopt("-c or --clic", "write output files needed by clic"); 
    tellopt("-d <x> or --deterministic <x>",
            "Deterministic algorithm, ignoring states with Pr <= x"); 
    tellopt("-F <x> or --scaleFactor <x>", "set DE scale factor");
    tellopt("-h or --help", "print this message");
    tellopt("-p <x> or --ptsPerDim <x>", "number of DE points per free var");
    tellopt("-s <x> or --strategy <x>", "set DE strategy");
    tellopt("-S <g>@<r> or --stage <g>@<r>",
            "add stage with <g> generations and <r> simulation reps."
            " The \"@<r>\"\n"
            "      portion can be omitted if -d or --deterministic are used.");
    tellopt("--stateIn <filename>",
            "read initial state from new-style file. Option may be repeated.");
    tellopt("--stateOut <filename>",
            "write final state to file");
    tellopt("-T <x> or --tol <x>", "termination criterion");
    tellopt("-t <x> or --threads <x>", "number of threads (default is auto)");
    tellopt("-v or --verbose", "verbose output");
    tellopt("--version", "Print version and exit");
    tellopt("-x <x> or --crossover <x>", "set DE crossover probability");
    exit(1);
}

int main(int argc, char **argv) {

    // Install handler for keyboard interrupts.
    signal(SIGINT, sighandle);
    signal(SIGTERM, handleSIGTERM);
    signal(SIGUSR1, sighandle);

    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
        {"clic", no_argument, 0, 'c'},
        {"deterministic", required_argument, 0, 'd'},
        {"threads", required_argument, 0, 't'},
        {"crossover", required_argument, 0, 'x'},
        {"scaleFactor", required_argument, 0, 'F'},
        {"strategy", required_argument, 0, 's'},
        {"stage", required_argument, 0, 'S'},
        {"tol", required_argument, 0, 'T'},
        {"ptsPerDim", required_argument, 0, 'p'},
        {"singletons", no_argument, 0, '1'},
        {"stateIn", required_argument, 0, 'z'},
        {"stateOut", required_argument, 0, 'y'},
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, 0, 'v'},
        {"version", no_argument, 0, 'V'},
        {NULL, 0, NULL, 0}
    };

    hdr("legofit: estimate population history");

    int i, j;
    time_t currtime = time(NULL);
    unsigned long pid = (unsigned long) getpid();
    double lo_twoN = 1.0, hi_twoN = 1e7;    // twoN bounds
    double lo_t = 0.0, hi_t = 1e7;  // t bounds
    int nThreads = 0;           // total number of threads
    int doSing = 0;             // nonzero means use singleton site patterns
    int write_clic_pts = 0;
    int status, optndx;
    long simreps = 1000000;
    char lgofname[200] = { '\0' };
    char patfname[200] = { '\0' };
    char stateOutName[200] = { '\0' };
    NameList *stateInNames = NULL;
    FILE *stateOut = NULL;

    // DiffEv parameters
    double F = 0.3;
    double CR = 0.8;
    double ytol = -1.0;         // stop when yspread <= ytol
    double ytol_default_stochastic = 3e-5;
    double ytol_default_deterministic = 1e-10;
    int strategy = 2;
    int ptsPerDim = 10;
    int verbose = 0;
    int deterministic = 0;
    int empty_reps = 0;
    SimSched *simSched = SimSched_new();

    // Ignore IdSet objects with probabilities <= improbable.
    extern long double improbable;

    assert(SimSched_nStages(simSched) == 0);

#if defined(__DATE__) && defined(__TIME__)
    printf("# Program was compiled: %s %s\n", __DATE__, __TIME__);
#endif
    printf("# Program was run: %s\n", ctime(&currtime));

    printf("# cmd:");
    for(i = 0; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');
    fflush(stdout);

    rngseed = currtime ^ pid;

    // command line arguments
    for(;;) {
        char *end;
        i = getopt_long(argc, argv, "cd:T:t:F:p:s:S:a:vx:1h", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case 'c':
            write_clic_pts = 1;
            break;
        case 'd':
            deterministic = 1;
            improbable = strtold(optarg, &end);
            if(*end != '\0') {
                fprintf(stderr,"\nIllegal argument to -d or --deterministic.\n"
                        "Can't parse %s as a long double.\n\n", optarg);
                usage();
            }
            break;
        case ':':
        case '?':
            usage();
            break;
        case 't':
            nThreads = strtol(optarg, NULL, 10);
            if(nThreads <= 0) {
                fprintf(stderr,"%s:%d: number of threads must be positive\n",
                        __FILE__,__LINE__);
                usage();
            }
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
                fprintf(stderr, "%s:%d: buffer overflow reading arg %s\n",
                        __FILE__, __LINE__, optarg);
                exit(EXIT_FAILURE);
            }
            g = r = b;

            // Parse string like "123@456". g is for diffev generations
            // and points to "123". r is for simulation replicates and
            // points to "456". It is legal to omit the "@456" part,
            // but only if -d or --deterministic is also used.
            (void) strsep(&r, "@");
            if(strlen(g) == 0)
                usage();
            long stageGen = strtol(g, NULL, 10);

            long stageRep;
            if(r == NULL || strlen(r) == 0) {
                empty_reps = 1;
                stageRep = 0;
            }else{
                stageRep = strtol(r, NULL, 10);
                if(stageRep == 0)
                    empty_reps = 1;
            }
            simreps = stageRep;

            if(stageGen < 0 || stageRep < 0) {
                fprintf(stderr,"Argument to -S or --stage can't have"
                        " negative integers: \"%s\" is illegal.\n",
                        optarg);
                exit(EXIT_FAILURE);
            }
            SimSched_append(simSched, stageGen, stageRep);
        }
            break;
        case 'v':
            verbose = 1;
            break;
        case 'V':
            return 0;
        case 'T':
            ytol = strtod(optarg, 0);
            break;
        case 'x':
            CR = strtod(optarg, 0);
            break;
        case 'y':
            status =
                snprintf(stateOutName, sizeof(stateOutName), "%s", optarg);
            if(status >= sizeof(stateOutName)) {
                fprintf(stderr, "%s:%d: buffer overflow\n",
                        __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            if( access( stateOutName, F_OK ) != -1 ) {
                // file already exists
                fprintf(stderr, "Warning: overwriting file %s\n", stateOutName);
            }
            stateOut = fopen(stateOutName, "w");
            if(stateOut == NULL) {
                fprintf(stderr, "%s:%d: can't open \"%s\" for output.\n",
                        __FILE__, __LINE__, stateOutName);
                exit(EXIT_FAILURE);
            }
            break;
        case 'z':
            stateInNames = NameList_append(stateInNames, optarg);
            CHECKMEM(stateInNames);
            break;
        case '1':
            doSing = 1;
            break;
        case 'h':
            usage();
            break;
        default:
            fprintf(stderr, "Can't parse option %c\n", i);
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
    snprintf(patfname, sizeof(patfname), "%s", argv[optind + 1]);
    assert(patfname[0] != '\0');

    if(!deterministic && empty_reps) {
        fprintf(stderr,"\nAt least one -S or --stage argument specified"
                " zero simulation replicates. In\n"
                "other words, there was no \"@\" or there was no positive"
                " integer after the \"@\".\n"
                "This is only legal when -d or --deterministic is also"
                " used.\n");
        exit(EXIT_FAILURE);
    }

    if(deterministic) {
        Network_init(DETERMINISTIC);

        if(ytol < 0.0)
            ytol = ytol_default_deterministic;

        if(0 == SimSched_nStages(simSched)) {
            // Default simulation schedule.
            // 1000 generations of 1 replicate
            SimSched_append(simSched, 1000, 1);
        }else if(1 < SimSched_nStages(simSched)) {
            fprintf(stderr,"\nDon't use the -S or --stage argument"
                    " more than once when using -d or --deterministic\n");
            exit(EXIT_FAILURE);
        }
    }else{
        Network_init(STOCHASTIC);

        if(ytol < 0.0)
            ytol = ytol_default_stochastic;

        if(0 == SimSched_nStages(simSched)) {
            // Default simulation schedule.
            // Stage 1: 200 DE generations of 1000 simulation replicates
            // Stage 2: 100 generations of 10000 replicates
            SimSched_append(simSched, 200, 1000);
            SimSched_append(simSched, 100, 10000);
            SimSched_append(simSched, 1000, simreps);
        }
    }

    SimSched_print(simSched, stdout);

    Bounds bnd = {
        .lo_twoN = lo_twoN,
        .hi_twoN = hi_twoN,
        .lo_t = lo_t,
        .hi_t = hi_t
    };

    void *network = Network_new(lgofname, bnd);
    LblNdx lblndx = Network_getLblNdx(network);

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, rngseed);
    rngseed += 1;

    int dim = Network_nFree(network); // number of free parameters
    if(dim == 0) {
        fprintf(stderr, "Error@%s:%d: no free parameters\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    int npts = dim * ptsPerDim;

    char *parnames[dim];
    for(i=0; i<dim; ++i) {
        parnames[i] = strdup(Network_getNameFree(network, i));
        CHECKMEM(parnames[i]);
    }

    // DiffEv state array is a matrix with a row for each point
    // and a column for each parameter.
    State *state=NULL;
    if(stateInNames) {
        // read States from files
        state = State_readList(stateInNames, npts, dim,
                               (const char * const *) parnames);
        CHECKMEM(state);
        if(npts != State_npoints(state)) {
            fprintf(stderr, "Revising npts from %d to %d\n",
                    npts, State_npoints(state));
            npts = State_npoints(state);
        }
        // Copy 0'th point in State object into GPTree, replacing
        // values specified in .lgo file.
        double x[dim];
        State_getVector(state, 0, dim, x);
        if(Network_setParams(network, dim, x)) {
            fprintf(stderr, "%s:%d: free params violate constraints\n",
                    __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        
        if(!Network_feasible(network, 1)) {
            fflush(stdout);
            dostacktrace(__FILE__, __LINE__, stderr);
            fprintf(stderr, "%s:%s:%d: stateIn points not feasible\n", __FILE__,
                    __func__,__LINE__);
            Network_printParStore(network, stderr);
            exit(EXIT_FAILURE);
        }
  } else {
        // de novo State
        state = State_new(npts, dim);
        CHECKMEM(state);
        for(i = 0; i < dim; ++i)
            State_setName(state, i, Network_getNameFree(network, i));
        for(i = 0; i < npts; ++i) {
            double x[dim];
            Network_initStateVec(network, i, dim, x, rng);
            State_setVector(state, i, dim, x);
        }
    }
    assert(state);

    if(nThreads == 0)
        nThreads = ceil(0.75 * getNumCores());
    if(nThreads > npts)
        nThreads = npts;

    if(deterministic) {
        printf("# algorithm          : %s\n", "deterministic");
        printf("# ignoring probs <=  : %Lg\n", improbable);
    }else {
        printf("# algorithm          : %s\n", "stochastic");
    }
    printf("# DE strategy        : %d\n", strategy);
    printf("#    F               : %lg\n", F);
    printf("#    CR              : %lg\n", CR);
    printf("#    tolerance       : %lg\n", ytol);
    printf("# nthreads           : %d\n", nThreads);
    printf("# lgo input file     : %s\n", lgofname);
    printf("# site pat input file: %s\n", patfname);
    printf("# free parameters    : %d\n", dim);
    printf("# points in DE swarm : %d\n", npts);
    if(stateInNames) {
        printf("# input state file(s):");
        NameList_print(stateInNames, stdout);
        putchar('\n');
    }
    if(stateOut)
        printf("# output state file  : %s\n", stateOutName);
    printf("# %s singleton site patterns.\n",
           (doSing ? "Including" : "Excluding"));
#if COST==KL_COST
    printf("# cost function      : %s\n", "KL");
#elif COST==LNL_COST
    printf("# cost function      : %s\n", "negLnL");
#else
# error "Unknown cost function"
#endif

    // Construct name for output file to which we will write
    // points for use in estimating Hessian matrix.
    char ptsfname[200] = {0};
    if(write_clic_pts){
        char a[200], b[200];
        strcpy(a, basename(lgofname));
        strcpy(b, basename(patfname));
        char *chrptr = strrchr(a, '.');
        if(chrptr)
            *chrptr = '\0';
        chrptr = strrchr(b, '.');
        if(chrptr)
            *chrptr = '\0';
        int version = 1;
        while(1) {
            // increment version number until we find one that hasn't
            // been used.
            status=snprintf(ptsfname, sizeof ptsfname, "%s-%s-%d.pts",
                            a, b, version);
            if(status >= sizeof ptsfname)
                DIE("buffer overflow");
            if( access( ptsfname, F_OK ) == -1 ) {
                // file name doesn't exist, so use this name
                break;
            }else{
                // increment version number and try again
                ++version;
            }
        }
        printf("# pts output file    : %s\n", ptsfname);
    }

    // Observed site pattern frequencies
    BranchTab *rawObs = parseOpf(patfname, &lblndx);
    if(doSing) {
        if(!BranchTab_hasSingletons(rawObs)) {
            fprintf(stderr, "%s:%d: Command line includes singletons "
                    "(-1 or --singletons)\n"
                    "    but none are present in \"%s\".\n",
                    __FILE__, __LINE__, patfname);
            exit(EXIT_FAILURE);
        }
    } else {
        if(BranchTab_hasSingletons(rawObs)) {
            fprintf(stderr, "%s:%d: Command line excludes singletons "
                    "(neither -1 nor --singletons)\n"
                    "    but singletons are present in \"%s\".\n",
                    __FILE__, __LINE__, patfname);
            exit(EXIT_FAILURE);
        }
    }
    BranchTab *obs = BranchTab_dup(rawObs);
    BranchTab_normalize(obs);

    // parameters for cost function
    CostPar costPar = {
        .obs = obs,
        .network = network,
        .nThreads = nThreads,
        .doSing = doSing,
        .simSched = simSched
    };

    // Number of parameters in quadratic model used to estimate
    // Hessian matrix: 1 intercept
    // dim linear terms
    // dim squared terms
    // (dim*(dim-1))/2 cross product terms
    unsigned nQuadPar = 1 + 2*dim + (dim*(dim-1))/2;

    // Number of points to use in fitting quadratic model
    unsigned nQuadPts = 10*nQuadPar;

    // parameters for Differential Evolution
    DiffEvPar dep = {
        .dim = dim,
        .ptsPerDim = ptsPerDim,
        .refresh = 2,           // how often to print a line of output
        .strategy = strategy,
        .nthreads = nThreads,
        .verbose = verbose,
        .seed = ((unsigned long) time(NULL)) - 1ul,
        .F = F,
        .CR = CR,
        .jobData = &costPar,
        .JobData_dup = CostPar_dup,
        .JobData_free = CostPar_free,
        .objfun = costFun,
        .threadData = NULL,
        .ThreadState_new = ThreadState_new,
        .ThreadState_free = ThreadState_free,
        .state = state,
        .simSched = simSched,
        .ytol = ytol,
        .pb = PointBuff_new(dim, nQuadPts)
    };

    double estimate[dim];
    double cost, yspread;

    printf("Initial parameter values\n");
    Network_printParStore(network, stdout);

    // Flush just before diffev so output file will be as complete as
    // possible while diffev is running.
    fflush(stdout);

    DEStatus destat = diffev(dim, estimate, &cost, &yspread, dep, rng);

    // Get mean site pattern branch lengths
    if(Network_setParams(network, dim, estimate)) {
        fprintf(stderr, "%s:%d: free params violate constraints\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    BranchTab *bt = get_brlen(network, simreps, doSing, rng);
    //    BranchTab_print(bt, stdout);

    const char *whyDEstopped;
    switch(destat) {
    case ReachedGoal:
        whyDEstopped = "reached_goal";
        break;
    case FinishedIterations:
        whyDEstopped = "finished_iterations";
        break;
    case Interrupted:
        whyDEstopped = "was_interrupted";
        break;
    default:
        whyDEstopped = "stopped_for_an_unknown_reason";
    }

    printf("DiffEv %s. cost=%0.5le spread=%0.5le",
           whyDEstopped, cost, yspread);
#if COST==LNL_COST
    printf("  relspread=%e", yspread / cost);
#endif
    putchar('\n');

    printf("Fitted parameter values\n");
    Network_printParStoreFree(network, stdout);

    // Put site patterns and branch lengths into arrays.
    unsigned npat = BranchTab_size(bt);
    tipId_t pat[npat];
    long double brlen[npat];
    BranchTab_toArrays(bt, npat, pat, brlen);

    // Determine order for printing lines of output
    unsigned ord[npat];
    orderpat(npat, ord, pat);

    printf("#%14s %10s\n", "SitePat", "BranchLen");
    char buff[100];
    for(j = 0; j < npat; ++j) {
        char buff2[100];
        snprintf(buff2, sizeof(buff2), "%s",
                 patLbl(sizeof(buff), buff, pat[ord[j]], &lblndx));
        printf("%15s %10.7Lf\n", buff2, brlen[ord[j]]);
    }

    if(stateOut) {
        State_print(state, stateOut);
        fclose(stateOut);
    }

    if(write_clic_pts) {

        double S = BranchTab_sum(rawObs); // sum of site pattern counts
        double entropy = BranchTab_entropy(obs); // -sum p ln(p)
        if(nQuadPts != PointBuff_size(dep.pb)) {
            fprintf(stderr,"Warning@%s:%d: expecting %u points; got %u\n",
                    __FILE__,__LINE__, nQuadPts, PointBuff_size(dep.pb));
            nQuadPts = PointBuff_size(dep.pb);
        }

        assert( access( ptsfname, F_OK ) == -1 );

        FILE *qfp = fopen(ptsfname, "w");
        if(qfp == NULL) {
            fprintf(stderr,"%s:%d: can't write file %s: using stdout\n",
                    __FILE__,__LINE__,ptsfname);
            qfp = stdout;
        }

        // First line contains dimensions. Each row contains dim+1 values:
        // first lnL, then dim parameter values.
        fprintf(qfp, "%u %d\n", nQuadPts, dim);

        // header
        fprintf(qfp, "%s", "lnL");
        for(i=0; i < dim; ++i)
            fprintf(qfp, " %s", Network_getNameFree(network, i));
        putc('\n', qfp);

        double lnL;

        // print estimate first
#if COST==KL_COST
        lnL = -S*(cost + entropy);
#else
        lnL = -cost;
#endif
        fprintf(qfp, "%0.18lg", lnL);
        for(i=0; i < dim; ++i)
            fprintf(qfp, " %0.18lg", estimate[i]);
        putc('\n', qfp);

        // print contents of PointBuff, omitting any that exactly
        // equal the point estimate
        while(0 != PointBuff_size(dep.pb)) {
            double c, par[dim];
            c = PointBuff_pop(dep.pb, dim, par);

            // Is current point the estimate? If so, skip it.
            int is_estimate=1;
            if(c != cost)
                is_estimate=0;
            for(i=0; is_estimate && i < dim; ++i) {
                if(par[i] != estimate[i])
                    is_estimate = 0;
            }
            if(is_estimate)
                continue;

#if COST==KL_COST
            lnL = -S*(c + entropy); // Kullback-Leibler cost function
#else
            lnL = -c;               // negLnL cost function
#endif
            fprintf(qfp, "%0.18lg", lnL);
            for(i=0; i < dim; ++i)
                fprintf(qfp, " %0.18lg", par[i]);
            putc('\n', qfp);
        }
        if(qfp != stdout) {
            fclose(qfp);
            fprintf(stderr,"%d points written to file %s\n",
                    nQuadPts, ptsfname);
        }
    }

    PointBuff_free(dep.pb);
    BranchTab_free(bt);
    BranchTab_free(rawObs);
    BranchTab_free(obs);
    gsl_rng_free(rng);
    Network_sanityCheck(network, __FILE__, __LINE__);
    Network_free(network);
    SimSched_free(simSched);
    State_free(state);
    binom_free();
    MatCoal_freeExterns();
    for(i=0; i<dim; ++i)
        free(parnames[i]);

    fprintf(stderr, "legofit is finished\n");

    return 0;
}
