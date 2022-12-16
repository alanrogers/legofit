/**
@file simpat.c
@page simpat
@brief Tabulate site pattern frequencies from data produced by ms2sim
or msprime.

# simpat: tabulates site patterns

Simpat reads data generated by ms2sim or msprime, tabulates counts of
nucleotide site patterns, and write the result to standard output.

# Usage

    Usage: simpat [options] [inputfile]
      Reads from stdin by default. Options may include:
       -V or --version
          Print version and exit
       -h or --help
          Print this message

# Example

`simpat` parses a file generated using `ms2sim` or `msprime`. The
first few lines of input should look like this:

    npops = 3
    pop sampsize
    x 2
    y 1
    z 1
    0 0 1 0 0
    0 0 1 0 0
    0 1 0 0 0

Line 1 gives the number of populations. These need not correspond to
the populations of `ms` or `msprime`. For example, `ms` might define
samples from the same population at different points in time, and you
might want to treat these as distinct samples within `simpat`.

Line 2 is a header for the 3 lines that follow.

In lines 3-5, the 1st column is a label for a population, and the 2nd
column specifies the number of haploid samples from that
population. The number of lines in this section of the input should
equal `npops`, as specified on line 1. The order in which these
samples are listed should correspond to the order in which they occur
in the lines of data that follow.

All remaining lines have the same format. Each provides data for a
single nucleotide site, and each consists of the same number of
fields, separated by white space. The 1st field labels the current
chromosome (or replicate). Each remaining field is a haploid genotype,
with 0 representing the ancestral allele and 1 the derived
allele. These are given in the order specified in the previous section
of output (lines 3-5 in the example above). For example, our example
implies that columns 2-3 refer to population x, 4 refers to y, and 5
to z.

In the output, site pattern "x:y" refers to the pattern in which the
derived allele is present haploid samples from "x" and "y" but not on
those from other populations. Here is the output of a run with 4
populations, x, y, n, and d:

    # simpat version 1.67
    # Including singleton site patterns.
    # Number of site patterns: 14
    # Nucleotide sites: 7298790
    # Sites used: 7298790
    #       SitePat             E[count]
                  x      1401902.0000000
                  y      1239084.0000000
                  n       589582.0000000
                  d       598305.0000000
                x:y      1074236.0000000
                x:n        21136.0000000
                x:d        21693.0000000
                y:n        38330.0000000
                y:d        31022.0000000
                n:d      1408276.0000000
              x:y:n        45407.0000000
              x:y:d        44519.0000000
              x:n:d       324614.0000000
              y:n:d       460684.0000000

The left column lists the site patterns that occur in the data. The
right column gives the expected count of each site pattern. These are
not necessarily integers, because they represent averages over all
possible subsamples consisting of a single haploid genome from each
population. They are integers here only because this simulation
modeled a single haploid sample from each population.

@copyright Copyright (c) 2018, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "binary.h"
#include "misc.h"
#include "simreader.h"
#include "typedefs.h"
#include "error.h"
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct Stack Stack;

/// Treat a vector of tipId_t values as a push-down stack.
struct Stack {
    int         dim, nused;
    tipId_t    *buff;           // not locally owned
};

static void usage(void);
static Stack *Stack_new(int dim, tipId_t buff[dim]);
static void Stack_free(Stack * stk);
static void Stack_push(Stack * self, tipId_t x);
static void generatePatterns(int bit, int npops, Stack * stk, tipId_t pat,
                             int doSing);

const char *useMsg =
    "\nUsage: simpat [options] [inputfile]\n"
    "  Reads from stdin by default. Options may include:\n";

/// Print usage message and die.
static void usage(void) {
    fputs(useMsg, stderr);
    tellopt("-V or --version", "Print version and exit");
    tellopt("-h or --help", "Print this message");
    exit(1);
}

/// This stack is local to this file. It provides a bounds-controlled
/// interface to an external array, which is passed as an argument, buff,
/// to Stack_new.
static Stack *Stack_new(int dim, tipId_t buff[dim]) {
    Stack      *self = malloc(sizeof(Stack));
    CHECKMEM(self);
    self->dim = dim;
    self->buff = buff;
    self->nused = 0;
    return self;
}

/// Frees the stack but not the underlying buffer.
static void Stack_free(Stack * stk) {
    free(stk);
}

/// Add an entry to the tail of the stack, checking bounds.
static void Stack_push(Stack * self, tipId_t x) {
    if(self->nused == self->dim) {
        fprintf(stderr, "%s:%s:%d ERR: buffer overflow\n",
                __FILE__, __func__, __LINE__);
        exit(EXIT_FAILURE);
    }
    self->buff[self->nused++] = x;
}

/// Call as generatePatterns(0, npops, stk, 0); Recursive function,
/// which generates all legal site patterns and pushes them onto a
/// stack.
static void
generatePatterns(int bit, int npops, Stack * stk, tipId_t pat, int doSing) {
    if(npops >= 8*sizeof(tipId_t)) {
        fprintf(stderr,"%s:%s:%d: %d is too many populations: max is %lu\n",
                __FILE__,__func__,__LINE__,
                npops, 8*sizeof(tipId_t) - 1);
        exit(EXIT_FAILURE);
    }
    const tipId_t unity = 1;
    if(bit == npops) {
        // Recursion stops here. If current pattern is
        // legal, then push it onto the stack. Then return.

        // Exclude patterns with all bits on, or all bits off.
        if(pat == 0 || pat == (unity << npops) - 1ULL)
            return;
        // Exclude singleton patterns unless "doSing" is true.
        if(!doSing && isPow2(pat))
            return;
        Stack_push(stk, pat);
        return;
    }
    tipId_t     on = 1UL << bit;
    generatePatterns(bit + 1, npops, stk, pat | on, doSing);    // curr bit on
    generatePatterns(bit + 1, npops, stk, pat, doSing); // curr bit off
}

int main(int argc, char **argv) {
    int         i, j, status, optndx, done;
    int         doSing = 1;     // nonzero means use singleton site patterns
    char        errbuff[1024] = { '\0' };
    FILE       *ifp = stdin;

    static struct option myopts[] = {
        // {char *name, int has_arg, int *flag, int val}
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'V'},
        {NULL, 0, NULL, 0}
    };

    // command line arguments
    for(;;) {
        i = getopt_long(argc, argv, "hV", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case '?':
            usage();
            break;
        case 'h':
            usage();
            break;
        case 'V':
            printf("simpat version %s\n", GIT_VERSION);
            return 0;
        default:
            usage();
        }
    }

    // There can be one additional argument, the name of the input
    // file.
    switch(argc - optind) {
    case 0:
        ifp = stdin;
        break;
    case 1:
        {
            char *fname = argv[optind];
            ifp = fopen(fname, "r");
            if(ifp==NULL) {
                fprintf(stderr,"%s:%d: can't open %s for input.\n",
                        __FILE__,__LINE__, fname);
                exit(EXIT_FAILURE);
            }
        }
        break;
    default:
        fprintf(stderr,"Only one input file is allowed.\n");
        usage();
    }

    SimReader *r=SimReader_new(ifp);
    if(r == NULL) {
        fprintf(stderr,"%s:%d: Can't read input\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    int n = SimReader_sampleDim(r);

    printf("# simpat version %s\n", GIT_VERSION);

    unsigned long npat = (1UL << n) - 2UL;  // number of site patterns
    if(!doSing)
        npat -= n;
    printf("# %s singleton site patterns.\n",
           (doSing ? "Including" : "Excluding"));
    printf("# Number of site patterns: %lu\n", npat);
    tipId_t     pat[npat];
    long double patCount[npat];
    int         lblsize = 100;
    char        lblbuff[lblsize];
    memset(patCount, 0, sizeof(patCount));

    {
        // Stack is a interface to array "pat".
        Stack      *stk = Stack_new(npat, pat);

        // Put site patterns into array "pat".
        generatePatterns(0, n, stk, 0, doSing);

        Stack_free(stk);
    }

    // Sort site patterns. Major sort is by number of "on" bits,
    // so that singleton patterns come first, then doubletons, etc.
    // Secondary sort is by order in which labels are listed
    // on the command line.
    qsort(pat, (size_t) npat, sizeof(pat[0]), compare_tipId);
    fflush(stdout);

    unsigned long nsites = 0, nbadaa = 0, nfixed = 0;
    long snpndx = -1;

    // Read data
    fprintf(stderr, "Doing %s pass through data to tabulate patterns..\n",
            "single");
    done=0;
    while(!done) {
        status = SimReader_next(r);
        switch(status) {
        case 0:
            ++nsites;
            break;
        case EOF:
            done=1;
            continue;
        case ALLELE_MISMATCH:
        case NO_ANCESTRAL_ALLELE:
            ++nbadaa;
            ++nsites;
            continue;
        default:
            // something wrong
            mystrerror_r(status, errbuff, sizeof errbuff);
            fprintf(stderr,"%s:%d: input error (%s)\n",
                    __FILE__,__LINE__, errbuff);
            exit(EXIT_FAILURE);
        }

        ++snpndx;

        // p and q are frequencies of derived and ancestral alleles
        double p[n], q[n];
        for(j = 0; j < n; ++j) {
            p[j] = SimReader_daf(r,j); // derived allele freq
            q[j] = 1.0 - p[j];
        }

        // Contribution of current snp to each site pattern.  Inner
        // loop considers each bit in current pattern.  If that bit is
        // on, multiply z by the derived allele frequency, p. If
        // that bit is off, multiply by q=1-p. In the end, z is Prod
        // p[j]^bit[j] * q[j]^(1-bit[j]) where bit[j] is the value (0
        // or 1) of the j'th bit.
        for(i = 0; i < npat; ++i) {
            tipId_t     pattern = pat[i];
            double      z = 1.0;
            for(j = 0; j < n; ++j) {
                if(pattern & 1u)
                    z *= p[j];
                else
                    z *= q[j];
                pattern >>= 1u;
            }
            if(!isfinite(z)) {
                fprintf(stderr, "%s:%d nonfinite z=%lf\n",
                        __FILE__, __LINE__, z);
                fprintf(stderr, "   pattern=%d\n", pat[i]);
                for(j = 0; j < n; ++j)
                    fprintf(stderr, "   %d: p=%lf q=%lf\n", j, p[j], q[j]);
            }
            assert(0 == (pattern & 1));
            patCount[i] += z;
        }
    }
    printf("# Nucleotide sites: %lu\n", nsites);
    if(nbadaa)
        printf("# Disagreements about alleles: %lu\n", nbadaa);
    if(nfixed)
        printf("# Monomorphic sites: %lu\n", nfixed);
    printf("# Sites used: %lu\n",
           nsites - nbadaa - nfixed);

    // patLbl needs an object of type LblNdx
    LblNdx lndx;
    LblNdx_init(&lndx);
    for(i=0; i < SimReader_sampleDim(r); ++i) {
        LblNdx_addSamples(&lndx,
                          1,
                          SimReader_lbl(r, i));
    }

    // print labels and site patterns
    printf("# %13s %20s", "SitePat", "E[count]");
    putchar('\n');
    for(i = 0; i < npat; ++i) {
        printf("%15s %20.7Lf\n",
               patLbl(lblsize, lblbuff, pat[i], &lndx), patCount[i]);
    }

    SimReader_free(r);
    if(ifp!=stdin)
        fclose(ifp);
    fprintf(stderr, "simpat is finished\n");
    return 0;
}
