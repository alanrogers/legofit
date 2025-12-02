/**
@file legosim.c
@page legosim
@brief Generate site patterns by coalescent simulation

# legosim: coalescent simulations within a network of populations

    usage: legosim [options] input_file_name
       where options may include:
       -b or --branch_length
          print branch lengths rather than probabilities
       -i <x> or --nItr <x>
          number of iterations in simulation
       -1 or --singletons
          Use singleton site patterns
       -d <x> or --deterministic <x>
          Deterministic algorithm, ignoring states with Pr <= x
       --network
          Print summary of population network and exit
       --plot <filename> or -p <filename>
          Create a dot-format file for plotting the network.
       -U <x>
          Mutations per generation per haploid genome.
       -h or --help
          print this message
       --version
          print version and exit
    Options -i, --nItr, and -U cannot be used with --deterministic.
    Option -U cannot be used with -b or --branch_length.

Here, "input_file" should be in @ref lgo ".lgo" format, which
describes the history of population size, subdivision, and gene
flow. The output looks like this:

    #######################################
    # legosim: site pattern probabilities #
    #   version 2.0.1-12-gf241700-dirty   #
    #######################################
    
    # Program was compiled: Nov 21 2020 10:15:39
    # Program was run: Sat Nov 21 10:17:57 2020
    
    # cmd: ./legosim -1 -d 0 input.lgo
    # input file                  : input.lgo
    # algorithm                   : deterministic
    # ignoring probs <=           : 0
    # including singleton site patterns.
    #       SitePat            Prob
                  x    0.1560242402
                  y    0.1481456533
                  n    0.4023561523
                x:y    0.2854101359
                x:n    0.0000926157
                y:n    0.0079712026

Here, the "SitePat" column labels site patterns. For example, site
pattern xy (denoted by "x:y" in this output) refers to nocleotide
sites at which the derived allele is present in single haploid samples
from X and Y but not in samples from other populations. The
probability of this site pattern, under the model of history specified
in file `input.lgo`, is given in the "Prob" column. When the `-d` or
`--deterministic` options are used, these probabilities are calculated
by a deterministic algorithm.

In the example above, I used `-d 0`, which tells legosim to ignore
only those states with probability 0. Thus, the answers will be
extremely accurate. When neither `-d` nor `--deterministic` are used,
the probabilities are estimated by coalescent simulation.

## Ignoring states of low probability

The deterministic algorithm sums across all possible histories of the
samples defined in the .lgo file. The number of histories increases
rapidly with sample size and with the number of migration events. For
larger models, it is not feasible to use `-d 0`. However, one can
still use the deterministic algorithm to obtain an approximate answer,
using an argument such as `-d 1e-6`. This tells legosim to use the
deterministic algorithm while ignoring states whose probability is
less than or equal to 1e-6. Nonzero arguments to `-d` will introduce
bias. To evaluate the magnitude of this bias, compare the results to
those obtained using the stochastic algorithm or (where feasible) the
deterministic one with `-d 0`.

When the deterministic algorithm is used, legosim prints a line of
output that looks like this:

    # Summed probability of ignored improbable events: 3.515617e-04

This "summed probability" can exceed 1, because it sums across several
probability distributions. Nonetheless, when it is small, that is an
indication that we have not ignored much of any of these distributions.
To rexpress this as a fraction of the of the total probability mass,
run legosim with argument `-d 1`. This ignores all events, so the
summed probability of ignored events will be an integer, which equals
the number of probability distributions involved. Divide the original
summed probability by this number to express it as a fraction of total
probability mass.

## Simulating site pattern counts

To simulate site pattern counts across an entire genome, use the `-U`
option, whose argument give the expected number of mutations per
generation per haploid genome. This argument should equal @f$\mu L@f$,
where @f$\mu@f$ is the mutation rate per nucleotide site per
generation, and @f$L@f$ is the number of nucleotide sites sequenced per
haploid genome, including monomorphic sites but excluding those that
fail quality control. For example, adding `-U 18` to the command above
led to the following output:

    #######################################
    # legosim: site pattern probabilities #
    #            version 1.89             #
    #######################################
    
    # Program was compiled: Oct 22 2020 22:30:28
    # Program was run: Fri Oct 23 12:35:53 2020
    
    # cmd: legosim -1 -i 10000 -U 18 input.lgo
    # input file                  : input.lgo
    # algorithm                   : stochastic
    # nreps                       : 10000
    # mutations per haploid genome: 18.000000
    # including singleton site patterns.
    #       SitePat           Count
                  x               4
                  y               1
                  n               8
                x:y               5
                x:n               0
                y:n               0

Now, the 2nd column is labeled "Count" rather than "Prob" and gives
the simulated count of each site pattern across the genome as a
whole. It is calculated by sampling from a Poisson distribution with
mean @f$\mu L b@f$, where @f$\mu@f$ and @f$L@f$ are as described
above, and @f$b@f$ is the average branch length. This Poisson model
treats nucleotide sites as independent, ignoring linkage
disequilibrium. The counts it provides are correct in expectation, but
their variances in repeated runs of the program are too small.

## Printing a summary of the network of nodes

Use the `--network` option to print a summary of the network of
nodes. As an example, let us use the following file, which is called
`toy.lgo`: 

    # Toy .lgo file
    time fixed    zero=0
    twoN fixed     one=1
    time free     Txyn=25000       # archaic-modern separation time
    time fixed    Tn=2000          # time of Neanderthal admixture
    time fixed    Txy=4000         # Africa-Eurasia separation time
    twoN free   twoNn=1e3          # archaic population size
    twoN free  twoNxy=1e4          # early modern population size
    mixFrac free  mN=0.02          # Neanderthal admixture into y
    segment x     t=zero   twoN=one    samples=1  # Africa
    segment y     t=zero   twoN=one    samples=1  # Eurasia
    segment n     t=Tn     twoN=twoNn  samples=1  # Neanderthal
    segment y2    t=Tn     twoN=twoNxy            # pre-mig eurasia
    segment xy    t=Txy    twoN=twoNxy            # early modern
    segment xyn   t=Txyn   twoN=twoNn             # ancestral
    mix    y  from y2 + mN * n      # y is a mixture of y2 and n
    derive x  from xy               # x is child of xy
    derive y2 from xy               # y2 is child of xy
    derive xy from xyn              # xy is child of xyn
    derive n  from xyn              # n is child of xyn

At the command line, type `legosim --network toy.lgo`. This produces
output that ends with:

    xyn: twoN=1000 ntrval=(25000.000000,inf)
      children: xy n
    |xy: twoN=10000 ntrval=(4000.000000,25000.000000)
    |  parents: xyn
    |  children: x y2
    ||x: twoN=1 ntrval=(0.000000,4000.000000)
    ||  parents: xy
    ||y2: twoN=10000 ntrval=(2000.000000,4000.000000)
    ||  parents: xy
    ||  children: y
    |||y: twoN=1 ntrval=(0.000000,2000.000000)
    |||  parents: y2 n
    |n: twoN=1000 ntrval=(2000.000000,25000.000000)
    |  parents: xyn
    |  children: y
    Samples:
      0     x
      1     y
      2     n


The first line of output tells us that xyn is the root segment, that
its population is `twoN=1000`, and that is spans an interval ranging
from 2500 generations ago to infinity. The second line tells us that
it has two "children", segments xy and n. The children and their
children are described in the lines that follow, preceded by "|"
characters to indicate the depth of nesting.

## Plotting the network

To make a graph of the network, use the `--plot <f.dot>` option, where
`f.dot` is the name of the file to be written as output. If that file
already exists, it will be overwritten. This output file is written in
the .dot format used by [GraphViz][GraphViz] to describe directed
graphs. I ran `legosim --plot toy.dot toy.lgo`, which produced
toy.lgo, an output file that looks like this

    digraph {
      ratio = 1.5;
      margin = 0;
      x [shape=square];
      y [shape=square];
      n [shape=square];
      { rank = same; x -> y [style = invis]};
      { rank = same; y2 n};
      xyn -> {xy n};
      xy -> {x y2};
      y2 -> y;
      n -> y;
    }

You can edit this file to tweak the resulting graph. In the first
line, `ratio = 1.5` sets the aspect ratio--the ratio of height to
width, and `margin = 0` says not to put white space around the graph.
The "shape=square" lines ensure that segments with data are plotted
with squares to distinguish them from segments without data. The line
`{ rank = same; x -> y [style = invis]}` says that x and y should have
the same "rank" (vertical position) and that x should be to the left
of right. Legosim emits a line like this so that all segments that end
at time zero have the same rank and are plotted in the order that they
appear in the .lgo file. The next line, `{ rank = same; y2 n}` says
that y2 and n are in the same rank but allows GraphViz to decide how
to arrange them within that rank. The remaining lines describe the
arrows to be drawn between segments.

To turn this input into a .png graphics file, type

    dot -Tpng toy.dot > toy.png

This produces the graph shown below.

![Graph of network in toy.lgo](../master/assets/toy.png)

I often find it useful to edit the .dot file, adding lines to place
segments in the same rank or to constrain their positions within a
rank. To do this, just follow the conventions illustrated in the .dot
file above.

[GraphViz]: https://graphviz.com

@copyright Copyright (c) 2025, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "branchtab.h"
#include "comb.h"
#include "gptree.h"
#include "lblndx.h"
#include "matcoal.h"
#include "network.h"
#include "parstore.h"
#include "patprob.h"
#include <assert.h>
#include <getopt.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

extern pthread_mutex_t seedLock;
extern unsigned long rngseed;
extern tipId_t union_all_samples;

void        usage(void);

void usage(void) {
    fprintf(stderr, "usage: legosim [options] input_file_name\n");
    fprintf(stderr, "   where options may include:\n");
    tellopt("-b or --branch_length", "print branch lengths rather than"
            " probabilities");
    tellopt("-i <x> or --nItr <x>", "number of iterations in simulation");
    tellopt("-1 or --singletons", "Use singleton site patterns");
    tellopt("-d <x> or --deterministic <x>",
            "Deterministic algorithm, ignoring states with Pr <= x");
    tellopt("--network", "Print summary of population network and exit");
    tellopt("--plot <f>", "Write dot-format plot commands to file f");
    tellopt("-U <x>", "Mutations per generation per haploid genome.");
    tellopt("-h or --help", "print this message");
    tellopt("--version", "print version and exit");
    fprintf(stderr,"Options -i, --nItr, and -U cannot be used with"
            " --deterministic.\n");
    fprintf(stderr,"Option -U cannot be used with"
            " -b or --branch_length.\n");
    exit(1);
}

int main(int argc, char **argv) {

    static struct option myopts[] =
        {
        /* {char *name, int has_arg, int *flag, int val} */
         {"branch_length", no_argument, 0, 'b'},
         {"deterministic", required_argument, 0, 'd'},
         {"network", no_argument, 0, 'n'},
         {"plot", required_argument, 0, 'p'},
         {"nItr", required_argument, 0, 'i'},
         {"mutations", required_argument, 0, 'U'},
         {"singletons", no_argument, 0, '1'},
         {"help", no_argument, 0, 'h'},
         {"version", no_argument, 0, 'V'},
         {NULL, 0, NULL, 0}
        };
    hdr("legosim: site pattern probabilities");

    int         i, j;
    int         doSing=0;  // nonzero => use singleton site patterns
    int         print_brlen=0; // nonzero => print branch lengths
    time_t      currtime = time(NULL);
    unsigned long pid = (unsigned long) getpid();
    double      lo_twoN = 1.0, hi_twoN = DBL_MAX;  // twoN bounds
    double      lo_t = 0.0, hi_t = DBL_MAX;        // t bounds
    double      U=0.0;          // mutations pre gen per haploid genome
    int         optndx;
    long        nreps = 100;
    int         deterministic = 0, stochastic = 0;
    int         show_network = 0;
    char        *plotfile = NULL;
    char        fname[200] = { '\0' };

    // Ignore IdSet objects with probabilities <= improbable.
    extern long double improbable, pr_ignored;

#if defined(__DATE__) && defined(__TIME__)
    printf("# Program was compiled: %s %s\n", __DATE__, __TIME__);
#endif
    printf("# Program was run: %s", ctime(&currtime));
    puts("#");
    printf("# cmd:");
    for(i = 0; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');
    char *dirname = getcwd(NULL, 0);
    printf("# curr dir: %s\n", dirname);
    free(dirname);

    // command line arguments
    for(;;) {
        char *end;
        i = getopt_long(argc, argv, "bd:np:i:U:1h", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            usage();
            break;
        case 'b':
            print_brlen = 1;
            break;
        case 'd':
            deterministic = 1;
            improbable = strtold(optarg, &end);
            if(*end != '\0') {
                fprintf(stderr,"Can't parse %s as a long double\n", optarg);
                exit(EXIT_FAILURE);
            }
            break;
        case 'n':
            show_network = 1;
            break;
        case 'p':
            plotfile = optarg;
            break;
        case 'e':
            break;
        case 'i':
            stochastic = 1;
            nreps = strtol(optarg, &end, 10);
            if(*end != '\0') {
                fprintf(stderr,"Can't parse %s as an integer\n", optarg);
                exit(EXIT_FAILURE);
            }
            break;
        case 'U':
            stochastic = 1;
            U = strtod(optarg,&end);
            if(*end != '\0') {
                fprintf(stderr,"Can't parse %s as a double\n", optarg);
                exit(EXIT_FAILURE);
            }
            break;
        case '1':
            doSing=1;
            break;
        case 'V':
            return 0;
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

    if(deterministic && stochastic) {
        fprintf(stderr,"\nOptions -d and --deterministic cannot be used"
                " with -i, --nItr, or -U.\n\n");
        usage();
    }

    if(print_brlen && U) {
        fprintf(stderr,"\nOptions -b and --branch_length cannot be used"
                " with -U.\n\n");
        usage();
    }

    if(deterministic)
        Network_init(DETERMINISTIC);
    else
        Network_init(STOCHASTIC);

    printf("# input file                  : %s\n", fname);
    if(show_network) {
        // no output
    }else if(deterministic) {
        printf("# algorithm                   : %s\n", "deterministic");
        printf("# ignoring probs <=           : %Lg\n", improbable);
        printf("# %s singleton site patterns.\n",
               (doSing ? "including" : "excluding"));
    }else {
        printf("# algorithm                   : %s\n", "stochastic");
        printf("# nreps                       : %lu\n", nreps);
        if(U)
            printf("# mutations per haploid genome: %lf\n", U);
        else
            printf("# not simulating mutations\n");
        printf("# %s singleton site patterns.\n",
               (doSing ? "including" : "excluding"));
    }

    Bounds bnd = {
            .lo_twoN = lo_twoN,
            .hi_twoN = hi_twoN,
            .lo_t = lo_t,
            .hi_t = hi_t
    };
    void *network = Network_new(fname, bnd);
    LblNdx lblndx = Network_getLblNdx(network);

    if(show_network) {
        Network_print(network, stdout);
        return 0;
    }

    if(plotfile) {
        printf(" writing plot commands to file %s\n", plotfile);
        FILE *fp = fopen(plotfile, "w");
        Network_plot(network, fp);
        fclose(fp);
        return 0;
    }

    int dim = Network_nFree(network);
    double x[dim];
    Network_getParams(network, dim, x);

    // No need to lock rngseed, because only 1 thread is running.
    rngseed = currtime^pid;
    gsl_rng  *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, rngseed);
    rngseed += 1;  // wraps to 0 at ULONG_MAX

    BranchTab *bt = get_brlen(network, nreps, doSing, 0.0, rng);

    if(!doSing)
        assert(!BranchTab_hasSingletons(bt));

    if(!U && !print_brlen)
        BranchTab_normalize(bt);

    if(deterministic)
        printf("# Summed probability of ignored improbable events: %Le\n",
               pr_ignored);

    // Put site patterns and branch lengths into arrays.
    unsigned npat = BranchTab_size(bt);
    tipId_t pat[npat];
    long double elen[npat];

    BranchTab_toArrays(bt, npat, pat, elen);
    //Network_printParStore(network, stdout);

    // Determine order for printing lines of output
    unsigned ord[npat];
    orderpat(npat, ord, pat);

    if(deterministic && improbable > 0 && print_brlen)
        fprintf(stderr, "\n# Warning: E[brlen] may be biased downward"
                " because the argument to -d\n"
                "# or --deterministic was > 0.\n#\n");

    if(U)
        printf("#%14s %15s\n", "SitePat", "Count");
    else if(print_brlen)
        printf("#%14s %15s\n", "SitePat", "E[brlen]");
    else
        printf("#%14s %15s\n", "SitePat", "Prob");
        
    char        buff[1000], buff2[1000];
    for(j = 0; j < npat; ++j) {
        if(!doSing && isPow2(pat[ord[j]]))
           continue;
        if(pat[ord[j]] == union_all_samples)
            continue;
        snprintf(buff2, sizeof(buff2), "%s",
                 patLbl(sizeof(buff), buff, pat[ord[j]], &lblndx));
        if(U) {
            unsigned mutations;
            mutations = gsl_ran_poisson(rng, U*elen[ord[j]]);
            printf("%15s %15u\n", buff2, mutations);
        }else
            printf("%15s %15.10Lf\n", buff2, elen[ord[j]]);
    }

    gsl_rng_free(rng);
    BranchTab_free(bt);
    Network_sanityCheck(network, __FILE__, __LINE__);
    Network_free(network);
    binom_free();
    MatCoal_freeExterns();

    return 0;
}
