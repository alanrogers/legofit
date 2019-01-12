/**
@file sitepat.c
@page sitepat
@brief Tabulate site pattern frequencies from .raf files.

# Sitepat: tabulates site patterns

Sitepat reads data in .raf format and tabulates counts of nucleotide
site patterns, writing the result to standard output. Optionally, it
also calculates a moving-blocks bootstrap, writing each bootstrap
replicate into a separate file.

# Usage

    Usage: sitepat [options] <x>=<in_1> <y>=<in_2> ... outgroup=<in_K>
       where <x> and <y> are arbitrary labels, and <in_i> are input
       files in raf format. Labels may not include the character ":".
       Final label must be "outgroup". Writes to standard output.
    
       If input file name ends with .gz, input is decompressed using
       gunzip.
     Maximum number of input files: 64 plus outgroup.
    
    Options may include:
       -f <name> or --bootfile <name>
          Bootstrap output file basename. Def: sitepat.boot.
       -r <x> or --bootreps <x>
          # of bootstrap replicates. Def: 0
       -b <x> or --blocksize <x>
          # of SNPs per block in moving-blocks bootstrap. Def: 0.
       -1 or --singletons
          Use singleton site patterns
       -m or --logMismatch
          Log REF mismatches to sitepat.log
       -A or --logAA
          Log sites with uncallable ancestral allele
       --version
          Print version and exit
       -h or --help
          Print this message

Reference genomes are often useful as outgroups, but they must first
be aligned to the target genome. These alignments are typically
distributed in axt format. To convert from axt format to raf format,
see @ref axt2raf "axt2raf".

# Example

Before running `sitepat`, use @ref raf "raf" to convert the input data
into raf format. To save space, you may want to compress the ".raf"
files using the external utility "gzip". The names of these compressed
files should end in ".gz", and `sitepat` will use an external utility
(gunzip) to decompress them. If `gunzip` has not been installed,
`sitepat` cannot read compressed input files. At present, gzipped
files cannot be used to generate bootstraps.

Let us assume you have generated raf files, and that directory ~/raf
contains a separate raf file for each population. We want to compare 4
populations, whose .raf files are `yri.raf`, `ceu.raf`, `altai.raf`,
and `denisova.raf`. Ancestral alleles are to be called using
`chimp.raf` as an outgroup. The following command will do this,
putting the results into `obs.txt`.

    sitepat x=raf/yri.raf \
           y=raf/ceu.raf \
           n=raf/altai.raf \
           d=raf/denisova.raf \
           outgroup=raf/chimp.raf > obs.txt

Here, "x", "y", "n", and "d" are labels that will be used to identify
site patterns in the output. For example, site pattern "x:y" refers to
the pattern in which the derived allele is present haploid samples
from "x" and "y" but not on those from other populations.  The order of
the command-line arguments determines the order in which labels are
sorted on output. Given the command line above, we would get a site
pattern labeled "x:y:d" rather than, say, "y:x:d".

The output looks like this:

    # Population labels:
    #    x = /home/rogers/raf/yri.raf
    #    y = /home/rogers/raf/ceu.raf
    #    n = /home/rogers/raf/altai.raf
    #    d = /home/rogers/raf/denisova.raf
    # Excluding singleton site patterns.
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

In the raf files used as input, chromosomes should appear in lexical
order. (See @ref sortorder "this link" for advice about lexical sorting.) 
Within each chromosome, nucleotides should appear in numerical  
order. There should be no duplicate (chromosome, position)
pairs. Otherwise, the program aborts with an error.

To generate a bootstrap, use the `--bootreps` option:

    sitepat --bootreps 50 \
           x=raf/yri.raf \
           y=raf/ceu.raf \
           n=raf/altai.raf \
           d=raf/denisova.raf \
           outgroup=raf/chimp.raf > obs.txt

This will generate not only the primary output file, `obs.txt`, but also
50 additional files, each representing a single bootstrap
replicate. The primary output file now has a bootstrap confidence
interval:

    # Population labels:
    #    x = /home/rogers/raf/yri.raf
    #    y = /home/rogers/raf/ceu.raf
    #    n = /home/rogers/raf/altai.raf
    #    d = /home/rogers/raf/denisova.raf
    # Excluding singleton site patterns.
    # Number of site patterns: 10
    # Tabulated 12327755 SNPs
    # bootstrap output file = sitepat.boot
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
interval. The bootstrap output files look like `sitepat.boot000`,
`sitepat.boot001`, and so on.

@copyright Copyright (c) 2016,2019 Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "binary.h"
#include "boot.h"
#include "rafreader.h"
#include "misc.h"
#include "strint.h"
#include "error.h"
#include "typedefs.h"
#include "version.h"
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXCHR 24               // maximum number of chromosomes

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
    "\nUsage: sitepat [options] <x>=<in_1> <y>=<in_2> ... outgroup=<in_K>\n"
    "   where <x> and <y> are arbitrary labels, and <in_i> are input\n"
    "   files in raf format. Labels may not include the character \":\".\n"
    "   Final label must be \"outgroup\". Writes to standard output.\n"
    "\n"
    "   If input file name ends with .gz, input is decompressed using\n"
    "   gunzip.\n";

/// Print usage message and die.
static void usage(void) {
    fputs(useMsg, stderr);
    fprintf(stderr, " Maximum number of input files: %lu plus outgroup.\n",
            8 * sizeof(bits_t));
    fputs("\nOptions may include:\n", stderr);
    tellopt("-f <name> or --bootfile <name>",
            "Bootstrap output file basename. Def: sitepat.boot.");
    tellopt("-r <x> or --bootreps <x>", "# of bootstrap replicates. Def: 0");
    tellopt("-b <x> or --blocksize <x>",
            "# of SNPs per block in moving-blocks bootstrap. Def: 0.");
    tellopt("-1 or --singletons", "Use singleton site patterns");
    tellopt("-m or --logMismatch", "Log REF mismatches to sitepat.log");
    tellopt("-A or --logAA", "Log sites with uncallable ancestral allele");
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

/// Add an entry to the stack, checking bounds.
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
    assert(npops < 8*sizeof(tipId_t));
    const tipId_t unity = 1;
    if(bit == npops) {
        // Recursion stops here. If current pattern is
        // legal, then push it onto the stack. Then return.

        // Exclude patterns with all bits on, or all bits off.
        if(pat == 0 || pat == (unity << npops) - unity)
            return;
        // Exclude singleton patterns unless "doSing" is true.
        if(!doSing && isPow2(pat))
            return;
        Stack_push(stk, pat);
        return;
    }
    tipId_t     on = unity << bit;
    generatePatterns(bit + 1, npops, stk, pat | on, doSing);    // curr bit on
    generatePatterns(bit + 1, npops, stk, pat, doSing); // curr bit off
}

int main(int argc, char **argv) {
    int         i, j, status, optndx, done;
    int         doSing = 0;     // nonzero means use singleton site patterns
    long        bootreps = 0;
    double      conf = 0.95;    // confidence level
    long        blocksize = 500;
    StrInt     *strint = StrInt_new();
    char        bootfname[FILENAMESIZE] = { '\0' };
    char        errbuff[100] = { '\0' };
    const char *logfname = "sitepat.log";
    int         logMismatch = 0, logAA = 0;
    FILE       *logfile = NULL;

    static struct option myopts[] = {
        // {char *name, int has_arg, int *flag, int val}
        {"bootfile", required_argument, 0, 'f'},
        {"bootreps", required_argument, 0, 'r'},
        {"blocksize", required_argument, 0, 'b'},
        {"singletons", no_argument, 0, '1'},
        {"logMismatch", no_argument, 0, 'm'},
        {"logAA", no_argument, 0, 'A'},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'V'},
        {NULL, 0, NULL, 0}
    };

    // command line arguments
    for(;;) {
        i = getopt_long(argc, argv, "b:c:f:hr:t:mAv1", myopts, &optndx);
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
        case 'V':
            printf("sitepat version %s\n", VERSION);
            return 0;
        case 'h':
            usage();
            break;
        case 'r':
            bootreps = strtol(optarg, NULL, 10);
            break;
        case '1':
            doSing = 1;
            break;
        case 'm':
            logMismatch = 1;
            break;
        case 'A':
            logAA = 1;
            break;
        default:
            usage();
        }
    }

    // remaining options: input files
    int         n = argc - optind;  // number of input files
    int         m = n-1;            // number excluding outgroup
    if(n == 0)
        usage();

    char       *poplbl[n];
    char       *fname[n];
    LblNdx      lndx;
    LblNdx_init(&lndx);
    RAFReader  *r[n];

    // Number of inputs can't exceed number of bits in an object of
    // type tipId_t.
    if(m > 8 * sizeof(tipId_t)) {
        fprintf(stderr, "Error: %d input files. Max is %lu.\n",
                n, 8*sizeof(tipId_t) + 1);
        usage();
    }
    // Parse remaining arguments, each of which should be of form
    // x=foo, where x is an arbitrary label and foo is the name of an
    // input file. Last label must be "outgroup".
    for(i = 0; i < n; ++i) {
        fname[i] = poplbl[i] = argv[i + optind];
        (void) strsep(fname + i, "=");
        if(fname[i] == NULL
           || poplbl[i] == NULL
           || strlen(poplbl[i]) == 0
           || strlen(fname[i]) == 0 || strchr(poplbl[i], ':') != NULL)
            usage();
        if(i < m)
            LblNdx_addSamples(&lndx, 1, poplbl[i]);
        r[i] = RAFReader_new(fname[i]);
    }
    if(0 != strcmp("outgroup", poplbl[n-1])) {
        fprintf(stderr,"%s:%d: last label is \"%s\""
                " instead of \"outgroup\".\n",
                __FILE__,__LINE__, poplbl[n-1]);
        usage();
    }

    if(logMismatch || logAA) {
        logfile = fopen(logfname, "w");
        if(logfile == NULL) {
            fprintf(stderr, "Can't write to file \"%s\".\n", logfname);
            exit(EXIT_FAILURE);
        }
    }

    // Default boot file name
    if(bootfname[0] == '\0') {
        const char *defName = "sitepat.boot";
        status = snprintf(bootfname, sizeof bootfname, "%s", defName);
        if(status >= sizeof bootfname) {
            fprintf(stderr, "%s:%d: ERR: Filename %s is too large."
                    " Max: %zu\n",
                    __FILE__, __LINE__, defName, sizeof(bootfname) - 1);
            exit(EXIT_FAILURE);
        }
    }

    printf("# sitepat version %s\n", VERSION);
    printf("# Population labels:\n");
    for(i = 0; i < n; ++i)
        printf("#  %s=%s\n", poplbl[i], fname[i]);

    // make sure labels are all different
    for(i = 1; i < n; ++i)
        for(j = 0; j < i; ++j)
            if(0 == strcmp(poplbl[i], poplbl[j])) {
                fprintf(stderr, "ERR: duplicate labels on command line.\n");
                fprintf(stderr, "     duplicated label: %s\n", poplbl[i]);
                exit(EXIT_FAILURE);
            }

    unsigned long npat = (1UL << m) - 2UL;  // number of site patterns
    if(!doSing)
        npat -= m;
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
        generatePatterns(0, m, stk, 0, doSing);
        Stack_free(stk);
    }

    // Sort site patterns. Major sort is by number of "on" bits,
    // so that doubleton patterns come first, then tripletons, ets.
    // Secondary sort is by order in which labels are listed
    // on the command line.
    qsort(pat, (size_t) npat, sizeof(pat[0]), compare_tipId);
    fflush(stdout);

    // Used by bootstrap
    Boot       *boot = NULL;
    int         nchr = 0;
    char        prev[RAFSTRSIZE], chr[RAFSTRSIZE] = { '\0' };
    long        nsnp[MAXCHR];
    memset(nsnp, 0, sizeof nsnp);

    // Read the data to get dimensions: number of chromosomes and
    // number of snps per chromosome. Then use these dimensions to
    // allocate a bootstrap object.
    if(bootreps > 0) {
        fprintf(stderr, "Doing 1st pass through data to get dimensions...\n");

        // First pass through data sets values of
        // nchr
        // nsnp[i] {i=0..nchr-1}
        done=0;
        while(!done) {
            status = RAFReader_multiNext(n, r);
            if(status==0)
                status = RAFReader_findDaf(n, r);
            switch(status) {
            case 0:
                break;
            case EOF:
                done=1;
                continue;
            case REF_MISMATCH:
            case MULTIPLE_ALT:
            case NO_ANCESTRAL_ALLELE:
                continue;
            default:
                // something wrong.
                mystrerror_r(status, errbuff, sizeof errbuff);
                fprintf(stderr,"%s:%d: input error (%s)\n",
                        __FILE__,__LINE__, errbuff);
                exit(EXIT_FAILURE);
            }

            assert(strlen(RAFReader_chr(r[0])) < sizeof prev);
            strcpy(prev, chr);
            strcpy(chr, RAFReader_chr(r[0]));
            int         diff = strcmp(prev, chr);
            if(diff != 0) {
                StrInt_insert(strint, chr, nchr);
                nsnp[nchr] = 1;
                ++nchr;
            } else
                ++nsnp[nchr - 1];
        }

        for(i = 0; i < n; ++i) {
            status = RAFReader_rewind(r[i]);
            if(status) {
                fprintf(stderr, "%s:%d: ERR: can't rewind input stream.\n",
                        __FILE__, __LINE__);
                fprintf(stderr, "  If --bootreps > 0, inputs must be"
                        " uncompressed files, not pipes.\n");
                exit(EXIT_FAILURE);
            }
        }

        // Allocate Boot structure
        gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set(rng, (unsigned long) time(NULL));
        boot = Boot_new(nchr, nsnp, bootreps, npat, blocksize, rng);
        gsl_rng_free(rng);
        CHECKMEM(boot);
    }

    unsigned long nsites = 0, nbadaa = 0, nbadref=0, nmultalt=0;
    long        snpndx = -1;

    // Iterate through raf files
    fprintf(stderr, "Doing %s pass through data to tabulate patterns..\n",
            bootreps > 0 ? "2nd" : "single");
    int         chrndx = -1, currChr = INT_MAX;
    RAFReader_clearChromosomes(n, r);
    done=0;
    while( !done ) {
        status = RAFReader_multiNext(n, r);
        if(status==0)
            status = RAFReader_findDaf(n, r);
        switch(status) {
        case 0:
            ++nsites;
            break;
        case EOF:
            done=1;
            continue;
        case REF_MISMATCH:
            ++nsites;
            ++nbadref;
            if(logMismatch) {
                fprintf(logfile,"REF mismatch:\n");
                RAFReader_printArray(n, r, logfile);
            }
            continue;
        case MULTIPLE_ALT:
            ++nsites;
            ++nmultalt;
            continue;
        case NO_ANCESTRAL_ALLELE:
            ++nsites;
            ++nbadaa;
            if(logAA) {
                fprintf(logfile,"Uncallable AA:\n");
                RAFReader_printArray(n, r, logfile);
            }
            continue;
        default:
            // something wrong.
            mystrerror_r(status, errbuff, sizeof errbuff);
            fprintf(stderr,"%s:%d: input error (%s)\n",
                    __FILE__,__LINE__, errbuff);
            exit(EXIT_FAILURE);
        }

        if(bootreps > 0) {
            // chrndx is index of current chromosome
            errno = 0;
            chrndx = StrInt_get(strint, RAFReader_chr(r[0]));
            if(errno) {
                fprintf(stderr,
                        "%s:%d: ERR: missing index for chromosome: %s\n",
                        __FILE__, __LINE__, RAFReader_chr(r[0]));
                exit(EXIT_FAILURE);
            }
            if(chrndx != currChr) {
                currChr = chrndx;
                snpndx = 0;
            } else
                ++snpndx;

#ifndef NDEBUG
            assert(snpndx < nsnp[chrndx]);
#endif
        }

        // p and q are frequencies of derived and ancestral alleles
        double      p[m], q[m];
        for(j = 0; j < m; ++j) {
            p[j] = RAFReader_daf(r[j]); // derived allele freq
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
            for(j = 0; j < m; ++j) {
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
                for(j = 0; j < m; ++j)
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
    printf("# Aligned sites                  : %lu\n", nsites);
    if(nbadref)
        printf("# Disagreements about ref allele : %lu\n", nbadref);
    if(nmultalt)
        printf("# Sites with multiple alt alleles: %lu\n", nmultalt);
    if(nbadaa)
        printf("# Undetermined ancestral allele  : %lu\n", nbadaa);
    printf("# Sites used                     : %lu\n",
           nsites - nbadaa - nbadref - nmultalt);

    // boottab[i][j] is the count of the j'th site pattern
    // in the i'th bootstrap replicate.
    double      bootvals[bootreps];
    double      boottab[bootreps][npat];
    memset(boottab, 0, sizeof boottab);

    if(bootreps > 0) {
        printf("# %s = %s\n", "bootstrap output file", bootfname);
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
            status = snprintf(buff, sizeof buff, "%s%d", bootfname, j);
            if(status >= sizeof buff)
                DIE("buffer overflow in snprintf");

            FILE       *fp = fopen(buff, "w");
            if(fp == NULL) {
                fprintf(stderr,"%s:%d: can't open \"%s\" for output.\n",
                        __FILE__,__LINE__,buff);
                exit(EXIT_FAILURE);
            }
            fprintf(fp, "# %13s %20s", "SitePat", "E[count]\n");
            for(i = 0; i < npat; ++i) {
                fprintf(fp, "%15s %20.7lf\n",
                        patLbl(lblsize, lblbuff, pat[i], &lndx),
                        boottab[j][i]);
            }
            fclose(fp);
        }
    }
    // print labels and binary representation of site patterns
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

    for(i = 0; i < n; ++i)
        RAFReader_free(r[i]);
    if(bootreps > 0)
        Boot_free(boot);
    StrInt_free(strint);
    if(logfile)
        fclose(logfile);
    fprintf(stderr, "sitepat is finished\n");
    return 0;
}
