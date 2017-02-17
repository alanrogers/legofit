/**
@file lego.c
@page lego
@brief Generate site patterns by coalescent simulation

# `lego`: coalescent simulations within a network of populations

    usage: lego [options] input_file
       where options may include:
       -i <x> or --nItr <x>
          number of iterations in simulation
       -1 or --singletons
          Use singleton site patterns
       -U <x>
          Mutations per generation per haploid genome.
       -h or --help
          print this message

Here, "input_file" should be in @ref lgo ".lgo" format, which
describes the history of population size, subdivision, and gene
flow. By default (i.e. if the `-U` option is not used), the output
looks like this:

    #########################################################
    # lego: generate site patterns by coalescent simulation #
    #########################################################

    # Program was compiled: Dec 25 2016 10:10:51
    # Program was run: Sun Dec 25 10:13:47 2016

    # cmd: ./lego -i 10000 input.lgo
    # nreps                       : 10000
    # input file                  : input.lgo
    # not simulating mutations
    # excluding singleton site patterns.
    #       SitePat E[BranchLength]
                x:y      39.7970280
                x:n      38.9656878
                y:n      40.8560014

Here, the "SitePat" column labels site patterns. For example, site
pattern xy (denoted by "x:y" in this output) refers to nocleotide
sites at which the derived allele is present in single haploid samples
from X and Y but not in samples from other populations. This site
pattern arises when a mutation strikes a branch that is ancestral only
to the samples from X and Y. The average length of this branch in
generations appears under "E[BranchLength]".

To simulate site pattern counts across an entire genome, use the `-U`
option, whose argument give the expected number of mutations per
generation per haploid genome. This argument should equal @f$\mu L@f$,
where @f$\mu@f$ is the mutation rate per nucleotide site per
generation, and @f$L@f$ is the number of nucleotide sites sequenced per
haploid genome, including monomorphic sites but excluding those that
fail quality control. For example, adding `-U 18` to the command above
led to the following output:

    #########################################################
    # lego: generate site patterns by coalescent simulation #
    #########################################################

    # Program was compiled: Dec 25 2016 10:10:51
    # Program was run: Mon Dec 26 09:31:59 2016

    # cmd: ./lego -i 10000 -U 18 input.lgo
    # nreps                       : 10000
    # input file                  : input.lgo
    # mutations per haploid genome: 18.000000
    # excluding singleton site patterns.
    #       SitePat           Count
                x:y             770
                x:n             713
                y:n             677

Now, the 2nd column is labeled "Count" rather than "E[BranchLength]"
and gives the simulated count of each site pattern across the genome
as a whole. It is calculated by sampling from a Poisson distribution
with mean @f$\mu L b@f$, where @f$\mu@f$ and @f$L@f$ are as described
above, and @f$b@f$ is the average branch length as reported without
the `-U` option. This Poisson model treats nucleotide sites as
independent, ignoring linkage disequilibrium. The counts it provides
are correct in expectation, but their variances in repeated runs of
the program are probably too small.

@copyright Copyright (c) 2015, 2016, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "gptree.h"
#include "patprob.h"
#include "parstore.h"
#include "lblndx.h"
#include "branchtab.h"
#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <limits.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

extern pthread_mutex_t seedLock;
extern unsigned long rngseed;

void        usage(void);

void usage(void) {
    fprintf(stderr, "usage: lego [options] input_file_name\n");
    fprintf(stderr, "   where options may include:\n");
    tellopt("-i <x> or --nItr <x>", "number of iterations in simulation");
	tellopt("-1 or --singletons", "Use singleton site patterns");
    tellopt("-U <x>", "Mutations per generation per haploid genome.");
    tellopt("-h or --help", "print this message");
    exit(1);
}

int main(int argc, char **argv) {

    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
        {"nItr", required_argument, 0, 'i'},
        {"mutations", required_argument, 0, 'U'},
        {"singletons", no_argument, 0, '1'},
        {"help", no_argument, 0, 'h'},
        {NULL, 0, NULL, 0}
    };

    printf("#########################################################\n"
           "# lego: generate site patterns by coalescent simulation #\n"
           "#########################################################\n");
    putchar('\n');

    int         i, j;
    int         doSing=0;  // nonzero => use singleton site patterns
    time_t      currtime = time(NULL);
	unsigned long pid = (unsigned long) getpid();
    double      lo_twoN = 0.0, hi_twoN = 1e6;  // twoN bounds
    double      lo_t = 0.0, hi_t = 1e6;        // t bounds
    double      U=0.0;          // mutations pre gen per haploid genome
    int         optndx;
    long        nreps = 100;
    char        fname[200] = { '\0' };
#if defined(__DATE__) && defined(__TIME__)
    printf("# Program was compiled: %s %s\n", __DATE__, __TIME__);
#endif
    printf("# Program was run: %s\n", ctime(&currtime));

    printf("# cmd:");
    for(i = 0; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');

    // command line arguments
    for(;;) {
        i = getopt_long(argc, argv, "i:t:U:1h", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            usage();
            break;
        case 'i':
            nreps = strtol(optarg, 0, 10);
            break;
        case 'U':
            U = strtod(optarg,NULL);
            break;
        case '1':
            doSing=1;
            break;
        case 'h':
        default:
            usage();
        }
    }

    // remaining option gives file name
    switch (argc - optind) {
    case 0:
        fprintf(stderr, "Command line must specify input file\n");
        usage();
        break;
    case 1:
        snprintf(fname, sizeof(fname), "%s", argv[optind]);
        break;
    default:
        fprintf(stderr, "Only one input file is allowed\n");
        usage();
    }
    assert(fname[0] != '\0');

    printf("# nreps                       : %lu\n", nreps);
    printf("# input file                  : %s\n", fname);
    if(U)
        printf("# mutations per haploid genome: %lf\n", U);
    else
        printf("# not simulating mutations\n");

    printf("# %s singleton site patterns.\n",
           (doSing ? "including" : "excluding"));

    Bounds bnd = {
            .lo_twoN = lo_twoN,
            .hi_twoN = hi_twoN,
            .lo_t = lo_t,
            .hi_t = hi_t
    };
    GPTree *gptree = GPTree_new(fname, bnd);
	LblNdx lblndx = GPTree_getLblNdx(gptree);

    int dim = GPTree_nFree(gptree);
    double x[dim];
    GPTree_getParams(gptree, dim, x);

    // No need to lock rngseed, because only 1 thread is running.
	rngseed = currtime^pid;
    gsl_rng  *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, rngseed);
	rngseed = (rngseed == ULONG_MAX ? 0 : rngseed+1);

    BranchTab *bt = patprob(gptree, nreps, doSing, rng);
    BranchTab_divideBy(bt, (double) nreps);
    //BranchTab_print(bt, stdout);

    // Put site patterns and branch lengths into arrays.
    unsigned npat = BranchTab_size(bt);
    tipId_t pat[npat];
    double prob[npat];
    double sqr[npat];
    BranchTab_toArrays(bt, npat, pat, prob, sqr);
	//GPTree_printParStore(gptree, stdout);

    // Determine order for printing lines of output
    unsigned ord[npat];
    orderpat(npat, ord, pat);

    if(U)
        printf("#%14s %15s\n", "SitePat", "Count");
    else
        printf("#%14s %15s\n", "SitePat", "E[BranchLength]");
    char        buff[100];
    for(j = 0; j < npat; ++j) {
        char        buff2[100];
        snprintf(buff2, sizeof(buff2), "%s",
                 patLbl(sizeof(buff), buff, pat[ord[j]], &lblndx));
        if(U) {
            unsigned mutations;
            mutations = gsl_ran_poisson(rng, U*prob[ord[j]]);
            printf("%15s %15u\n", buff2, mutations);
        }else
            printf("%15s %15.7lf\n", buff2, prob[ord[j]]);
    }

    GPTree_sanityCheck(gptree, __FILE__, __LINE__);
    GPTree_free(gptree);

    return 0;
}
