/**
@file scrmpat.c
@page scrmpat
@brief Tabulate site pattern frequencies from .daf files.

# Scrmpat: tabulates site patterns

Scrmpat reads data generated by scrm (with option -transpose-segsites)
tabulates counts of nucleotide site patterns, writing the result to
standard output. Optionally, it also calculates a moving-blocks
bootstrap, writing each bootstrap replicate into a separate file.

# Usage

    Usage: scrmpat [options] <x> <y> ...
      where <x>, <y>, etc. are arbitrary labels, whose number and order
      must agree with that of the populations specified in the scrm command
      line (using scrm arguments -I and -eI). Labels may not include the
      character ":". Writes to standard output. Max number of input files: 32.

    Options may include:
       --infile <name>
          Input file name. Def: standard input
       --bootfile <name>
          Bootstrap output file basename. Def: scrmpatX.bsr
       -r <x> or --bootreps <x>
          # of bootstrap replicates. Def: 0
       -b <x> or --blocksize <x>
          # of SNPs per block in moving-blocks bootstrap. Def: 0.
       -F or --logFixed
          log fixed sites to scrmpat.log
       -a or --logAll
          log all sites to scrmpat.log
       --version
          Print version and exit
       -h or --help
          Print this message

# Example

`scrmpat` parses a file generated using `scrm`. The `scrm` command
should include the option `-transpose-segsites`. Let us assume you
have done this, that file `foo.scrm` contains the output simulated
by `scrm`, and that these simulated data included genotypes referring
to four populations, labeled "x", "y", "n", and "d". The `scrmpat`
command would look like this:

    scrmpat --infile foo.scrm x y n d

`scrmpat`'s notion of a "population" differs from that of `scrm`, in
that `scrmpat` treats samples of different ages as separate
populations, even if they reside in the same population on the `scrm`
command line. For example, consider the following `scrm` command line:

    scrm 3 -I 2 1 1 -eI 0.5 0 1

This specifies three haploid samples distributed across two
populations. The `-I` argument says that each population has a
sample at time 0. The `-eI` argument says that, in addition,
population 2 has a sample at time 0.5. All three samples would be
treated as separate populations by `scrmpat`. Thus, the `scrmpat`
command line should list three labels, as in "scrmpat x y z".

In the output, site pattern "x:y" refers to
the pattern in which the derived allele is present haploid samples
from "x" and "y" but not on those from other populations. The order of
the command-line arguments determines the order in which labels are
sorted on output. Given the command line above, we would get a site
pattern labeled "x:y:d" rather than, say, "y:x:d".

The output looks like this:

    # scrmpat version 1.3
    # Population labels: x y n d
    # Number of site patterns: 10
    # Tabulated 12327755 SNPs
    #       SitePat             E[count]
                x:y       340952.4592501
                x:n        46874.1307236
                x:d        46034.4670204
                y:n        55137.4236715
                y:d        43535.5248078
                n:d       231953.3372578
              x:y:n        91646.1277991
              x:y:d        88476.9619569
              x:n:d        96676.3877423
              y:n:d       100311.4411513

The left column lists the site patterns that occur in the data. The
right column gives the expected count of each site pattern. These are
not integers, because they represent averages over all possible
subsamples consisting of a single haploid genome from each
population.

To generate a bootstrap, use the `--bootreps` option:

    scrmpat --bootreps 50 --infile foo.scrm x y n d > obs.txt

This will generate not only the primary output file, `obs.txt`, but also
50 additional files, each representing a single bootstrap
replicate. The primary output file now has a bootstrap confidence
interval:

    # Including singleton site patterns.
    # Number of site patterns: 10
    # Tabulated 12327755 SNPs
    # bootstrap output file = scrmpatX.bsr
    # confidence level = 95%
    #       SitePat             E[count]          loBnd          hiBnd
                x:y       340952.4592501 338825.6604586 342406.6670816
                x:n        46874.1307236  46361.5798377  47438.1857029
                x:d        46034.4670204  45605.6588012  46631.6434277
                y:n        55137.4236715  54650.0763578  55783.7051253
                y:d        43535.5248078  43110.5119922  44234.0919024
                n:d       231953.3372578 229495.3741057 234173.6878092
              x:y:n        91646.1277991  90494.0219749  92873.4443706
              x:y:d        88476.9619569  87137.1867967  89585.8431419
              x:n:d        96676.3877423  95935.5184294  97417.6241185
              y:n:d       100311.4411513  99292.9839140 101163.3457462

Here, `loBnd` and `hiBnd` are the limits of a 95% confidence
interval. The bootstrap output files look like `scrmpat0.bsr`,
`scrmpat1.bsr`, and so on.

@copyright Copyright (c) 2018, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "binary.h"
#include "boot.h"
#include "misc.h"
#include "scrmreader.h"
#include "typedefs.h"
#include "version.h"
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
    "\nUsage: scrmpat [options] <x> <y> ...\n"
    "  where <x>, <y>, etc. are arbitrary labels, whose number and order\n"
    "  must agree with that of the populations specified in the scrm command\n"
    "  line (using scrm arguments -I and -eI). Labels may not include the\n"
    "  character \":\". Writes to standard output.";

/// Print usage message and die.
static void usage(void) {
    fputs(useMsg, stderr);
    fputs("\nOptions may include:\n", stderr);
    tellopt("--infile <name>",
            "Input file name. Def: standard input");
    tellopt("--bootfile <name>",
            "Bootstrap output file basename. Def: scrmpatX.bsr");
    tellopt("-r <x> or --bootreps <x>", "# of bootstrap replicates. Def: 0");
    tellopt("-b <x> or --blocksize <x>",
            "# of SNPs per block in moving-blocks bootstrap. Def: 0.");
    tellopt("-F or --logFixed", "log fixed sites to scrmpat.log");
    tellopt("-a or --logAll", "log all sites to scrmpat.log");
    tellopt("--version", "Print version and exit");
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
    assert(sizeof(tipId_t) < sizeof(unsigned long long));
    if(bit == npops) {
        // Recursion stops here. If current pattern is
        // legal, then push it onto the stack. Then return.

        // Exclude patterns with all bits on, or all bits off.
        if(pat == 0 || pat == (1ULL << npops) - 1ULL)
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
    long        bootreps = 0;
    double      conf = 0.95;    // confidence level
    long        blocksize = 500;
    char        bootfname[FILENAMESIZE] = { '\0' };
    char        errbuff[100] = { '\0' };
    const char *logfname = "scrmpat.log";
    int         logFixed = 0, logAll = 0;
    FILE       *logfile = NULL;
    FILE       *ifp = stdin;

    static struct option myopts[] = {
        // {char *name, int has_arg, int *flag, int val}
        {"infile", required_argument, 0, 'i'},
        {"bootfile", required_argument, 0, 'f'},
        {"bootreps", required_argument, 0, 'r'},
        {"blocksize", required_argument, 0, 'b'},
        {"logFixed", no_argument, 0, 'F'},
        {"logAll", no_argument, 0, 'a'},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'V'},
        {NULL, 0, NULL, 0}
    };

    // command line arguments
    for(;;) {
        i = getopt_long(argc, argv, "ab:c:hi:r:t:Fv", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            usage();
            break;
        case 'b':
            blocksize = strtod(optarg, NULL);
            if(blocksize <= 0) {
                fprintf(stderr,
                        "%s:%d: bad argument to -b or --blocksize: \"%s\"\n",
                        __FILE__, __LINE__, optarg);
                usage();
            }
            break;
        case 'f':
            status = snprintf(bootfname, sizeof bootfname, "%s", optarg);
            if(status >= sizeof bootfname) {
                fprintf(stderr, "%s:%d: ERR: Filename %s is too large."
                        " Max: %zu\n",
                        __FILE__, __LINE__, optarg, sizeof(bootfname) - 1);
                exit(EXIT_FAILURE);
            }
            break;
        case 'h':
            usage();
            break;
        case 'i':
            ifp = fopen(optarg, "r");
            if(ifp==NULL) {
                fprintf(stderr,"%s:%d: can't open %s for input.\n",
                        __FILE__,__LINE__, optarg);
                exit(EXIT_FAILURE);
            }
            break;
        case 'V':
            printf("scrmpat version %s\n", VERSION);
            return 0;
        case 'r':
            bootreps = strtol(optarg, NULL, 10);
            break;
        case 'F':
            logFixed = 1;
            break;
        case 'a':
            logAll = 1;
            break;
        default:
            usage();
        }
    }

    // remaining options: population labels
    int         n = argc - optind;  // number of population labels
    if(n == 0)
        usage();

    char       *poplbl[n];
    LblNdx      lndx;
    LblNdx_init(&lndx);

    // Number of inputs can't exceed number of bits in an object of
    // type tipId_t.
    if(n > 8 * sizeof(tipId_t)) {
        fprintf(stderr, "Error: %d populations. Max is %lu.\n",
                n, 8 * sizeof(tipId_t));
        usage();
    }
    // Parse remaining arguments, each of which should be an arbitrary
    // label.
    for(i = 0; i < n; ++i) {
        poplbl[i] = argv[i + optind];
        if(poplbl[i] == NULL
           || strlen(poplbl[i]) == 0
           || strchr(poplbl[i], ':') != NULL)
            usage();
        LblNdx_addSamples(&lndx, 1, poplbl[i]);
    }

    if(logFixed || logAll) {
        logfile = fopen(logfname, "w");
        if(logfile == NULL) {
            fprintf(stderr, "Can't write to file \"%s\".\n", logfname);
            exit(EXIT_FAILURE);
        }
    }
    if(ifp==stdin && (bootreps>0 || bootfname[0] != '\0')) {
        fprintf(stderr, "%s:%d: Can't do bootstrap when input is stdin.\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    // Default boot file name
    if(bootfname[0] == '\0') {
        const char *defName = "scrmpat";
        status = snprintf(bootfname, sizeof bootfname, "%s", defName);
        if(status >= sizeof bootfname) {
            fprintf(stderr, "%s:%d: ERR: Filename %s is too large."
                    " Max: %zu\n",
                    __FILE__, __LINE__, defName, sizeof(bootfname) - 1);
            exit(EXIT_FAILURE);
        }
    }

    printf("# scrmpat version %s\n", VERSION);
    printf("# Population labels:");
    for(i = 0; i < n; ++i)
        printf(" %s", poplbl[i]);
    putchar('\n');

    // make sure labels are all different
    for(i = 1; i < n; ++i)
        for(j = 0; j < i; ++j)
            if(0 == strcmp(poplbl[i], poplbl[j])) {
                fprintf(stderr, "ERR: duplicate labels on command line.\n");
                fprintf(stderr, "     duplicated label: %s\n", poplbl[i]);
                exit(EXIT_FAILURE);
            }

    unsigned long npat = (1UL << n) - 2UL;  // number of site patterns
    if(!doSing)
        npat -= n;
    printf("# %s singleton site patterns.\n",
           (doSing ? "Including" : "Excluding"));
    printf("# Number of site patterns: %lu\n", npat);
    tipId_t     pat[npat];
    double      patCount[npat];
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

    // Used by bootstrap
    Boot       *boot = NULL;
    ScrmReader *r=ScrmReader_new(ifp);
    if(r == NULL) {
        fprintf(stderr,"%s:%d: Can't read scrm output\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    if(n != ScrmReader_sampleDim(r)) {
        fprintf(stderr,"%s:%d:"
                " Number (%d) of labels != dimension (%d) of sample array\n"
                "in scrm output.\n",
                __FILE__,__LINE__,n,ScrmReader_sampleDim(r));
        exit(EXIT_FAILURE);
    }

    // Read the data to get dimensions. nchr=1 by fiat: rather than
    // treating the independent replicates produced by scrm as
    // separate chromosomes, I concatenate them all into a single
    // chromosome and do a moving-blocks bootstrap on that. This
    // is because in simulations, chromosomes are often pretty short,
    // and the number of chromosomes is large. I'm afraid that a
    // moving-blocks bootstrap on short chromosomes would be poorly
    // behaved. So this code sets nchr=1, calculates the total number
    // of snps, and then uses these dimensions to allocate a bootstrap
    // object.
    int         nchr = 1;
    long        nsnp[1];
    nsnp[0] = 0;
    if(bootreps > 0) {
        fprintf(stderr, "Doing 1st pass through data to get dimensions...\n");

        // First pass through data sets values of
        // nchr, nsnp[0]
        done = 0;
        while(!done) {
            status = ScrmReader_next(r);
            switch(status) {
            case 0:
                break;
            case EOF:
                done=1;
                continue;
            case ALLELE_MISMATCH:
            case NO_ANCESTRAL_ALLELE:
                continue;
            default:
                // something wrong
                mystrerror_r(status, errbuff, sizeof errbuff);
                fprintf(stderr,"%s:%d: input error (%s)\n",
                        __FILE__,__LINE__, errbuff);
                exit(EXIT_FAILURE);
            }

            ++nsnp[0];
        }

        status = ScrmReader_rewind(r);
        if(status) {
            fprintf(stderr, "%s:%d: ERR: can't rewind input stream.\n",
                    __FILE__, __LINE__);
            fprintf(stderr, "  If --bootreps > 0, inputs must be"
                    " files, not pipes.\n");
            exit(EXIT_FAILURE);
        }

        // Allocate Boot structure
        gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set(rng, (unsigned long) time(NULL));
        boot = Boot_new(nchr, nsnp, bootreps, npat, blocksize, rng);
        gsl_rng_free(rng);
        CHECKMEM(boot);
    }

    unsigned long nsites = 0, nbadaa = 0, nfixed = 0;
    long snpndx = -1;
    int chrndx=0;

    // Read data
    fprintf(stderr, "Doing %s pass through data to tabulate patterns..\n",
            bootreps > 0 ? "2nd" : "single");
    done=0;
    while(!done) {
        status = ScrmReader_next(r);
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
            p[j] = ScrmReader_daf(r,j); // derived allele freq
            q[j] = 1.0 - p[j];
        }

        if(logAll) {
            fprintf(logfile, "%5u %10lu\n", ScrmReader_chr(r),
                    ScrmReader_nucpos(r));
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
            if(bootreps > 0) {
                assert(snpndx >= 0);
                assert(chrndx >= 0);
                Boot_add(boot, chrndx, snpndx, i, z);
            }

        }
#ifndef NDEBUG
        if(bootreps > 0)
            Boot_sanityCheck(boot, __FILE__, __LINE__);
#endif
    }
    printf("# Nucleotide sites: %lu\n", nsites);
    if(nbadaa)
        printf("# Disagreements about alleles: %lu\n", nbadaa);
    if(nfixed)
        printf("# Monomorphic sites: %lu\n", nfixed);
    printf("# Sites used: %lu\n",
           nsites - nbadaa - nfixed);

    // boottab[i][j] is the count of the j'th site pattern
    // in the i'th bootstrap replicate.
    double      bootvals[bootreps];
    double      boottab[bootreps][npat];
    memset(boottab, 0, sizeof boottab);

    if(bootreps > 0) {
        printf("# %s = %sX.bsr\n", "bootstrap output file", bootfname);
        printf("# %s = %4.2lf%%\n", "confidence level", 100 * conf);
#ifndef NDEBUG
        Boot_sanityCheck(boot, __FILE__, __LINE__);
#endif
        // put site pattern counts into matrix boottab.
        for(i = 0; i < bootreps; ++i)
            Boot_aggregate(boot, i, npat, boottab[i]);

        // write an output file for each bootstrap replicate
        for(j = 0; j < bootreps; ++j) {
            char        buff[FILENAMESIZE + 3];
            status = snprintf(buff, sizeof buff, "%s%d.bsr", bootfname, j);
            if(status >= sizeof buff)
                DIE("buffer overflow in snprintf");

            FILE       *fp = fopen(buff, "w");
            if(fp == NULL)
                DIE("bad fopen");
            fprintf(fp, "# %13s %20s", "SitePat", "E[count]\n");
            for(i = 0; i < npat; ++i) {
                fprintf(fp, "%15s %20.7lf\n",
                        patLbl(lblsize, lblbuff, pat[i], &lndx),
                        boottab[j][i]);
            }
            fclose(fp);
        }
    }
    // print labels and site patterns
    printf("# %13s %20s", "SitePat", "E[count]");
    if(bootreps > 0)
        printf(" %15s %15s", "loBnd", "hiBnd");
    putchar('\n');
    for(i = 0; i < npat; ++i) {
        printf("%15s %20.7lf",
               patLbl(lblsize, lblbuff, pat[i], &lndx), patCount[i]);
        if(bootreps > 0) {
            double      lowBnd, highBnd;
            for(j = 0; j < bootreps; ++j)
                bootvals[j] = boottab[j][i];
            confidenceBounds(&lowBnd, &highBnd, conf, bootreps, bootvals);
            printf(" %15.7lf %15.7lf", lowBnd, highBnd);
        }
        putchar('\n');
    }

    ScrmReader_free(r);
    if(bootreps > 0)
        Boot_free(boot);
    if(logfile)
        fclose(logfile);
    if(ifp!=stdin)
        fclose(ifp);
    fprintf(stderr, "scrmpat is finished\n");
    return 0;
}
