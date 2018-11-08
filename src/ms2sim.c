/**
 * @file ms2sim.c
 * @page ms2sim
 * @brief ms2sim translates output from Dick Hudson's `ms` program
 * into "sim" format.
 *
 * # ms2sim, a program that converts `ms` output into sim format
 *
 * This program is written in C, because the conversion involves
 * transposing a large matrix, and this is faster in C.  If you run it
 * like "ms2sim my_input_file_name", it will allocate internal arrays
 * using default values. If your `ms` output is large, these arrays will
 * be too small. The program will print an error message and abort. To
 * avoid this problem, run the (unix or linux) command "grep positions
 * ms.out | wc". This will print 3 numbers. The number of tokens must be
 * as large as 2nd; buffer size as large as the 3rd. You can set these
 * values using the `--maxTokens` and `--inBuffSize` arguments.
 *
 * For example, I used `ms` to produce a file called `ms.out`. To find out
 * how large to set the arrays, I executed
 *
 *     grep positions ms.out | wc
 *
 * which produced the output
 *
 *     1   89366 1161757
 *
 * This indicates that I need to accomodate input lines with 89366 tokens
 * and 1161757 characters. So I ran `ms2sim` like this:
 *
 *     ms2sim --maxTokens 89370 --inBuffSize 1161800 ms.out > ms.gtp
 *
 * This reads input from file `ms.out` and writes to file `ms.gtp`.
 * The first few lines of `ms.gtp` look like this:
 *
 *     # ms2sim was run Sat Feb  2 11:15:05 2013
 *     # cmd line     = ms2sim -t 89370 -b 1161800 ms.out
 *     # sim cmd      = ms 50 1 -t 20000 -r 20000 1000000
 *     # segsites     = 89365
 *     # nsequences   = 50
 *     # recomb rate per site = 1e-06
 *     # ploidy       = 1
 *     #   snp_id     nucpos         mappos alleles genotypes
 *              0          1 0.000000006802      01 0011000000...
 *              1          6 0.000000061255      01 0000000000...
 *
 * In this example, I've truncated the genotypes on the right in order
 * to fit them on the page.
 *
 * Here is the usage message:
 *
 *     usage: ms2sim [options] [input_file_name]
 *        where options may include:
 *        -t \<x\> or --maxTokens \<x\>
 *           maximum tokens per input line
 *        -b \<x\> or --inBuffSize \<x\>
 *           size of input buffer
 *        -R \<x\> or --recombination \<x\>
 *           rate for adjacent nucleotides
 *        -h     or --help
 *           print this message
 *       To figure out appropriate values for -t and -b, execute the
 *       following at the command line:
 *
 *          grep positions ms.out | wc
 *
 *       where "ms.out" is the ms output file. This command prints 3 numbers.
 *       Number of tokens must be as large as 2nd; buffer size as large as the
 *       3rd. If no input file is given, ms2sim reads stdin. ms2sim
 *       always writes to standard output.
 *
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "misc.h"
#include "tokenizer.h"
#include "uintqueue.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <getopt.h>
#include <assert.h>
#include <math.h>

typedef struct MSSample MSSample;

struct MSSample {
    double age;
    int pop;      // population
    int nsamp;    // number of haploid samples
    MSSample *next;
};

void      usage(void);
unsigned *countSamples(Tokenizer *tkz, int *sampleDim);
MSSample *MSSample_new(double age, int pop, int nsamp);
int       MSSample_compare(MSSample *left, MSSample *right);
MSSample *MSSample_insert(MSSample *node, MSSample *new);
MSSample *MSSample_collapse(MSSample *node);
int       MSSample_length(MSSample *self);

MSSample *MSSample_new(double age, int pop, int nsamp) {
    MSSample *self = malloc(sizeof(MSSample));
    CHECKMEM(self);
    memset(self, 0, sizeof(MSSample));
    return self;
}

void MSSample_free(MSSample *self) {
    if(self == NULL)
        return;
    MSSample_free(self->next);
    free(self);
}

// Sort by age.
int MSSample_compare(MSSample *left, MSSample *right) {
    if(left->age > right->age)
        return 1;
    if(left->age < right->age)
        return -1;
    // ages equal
    return 0;
}

int MSSample_length(MSSample *self) {
    if(self == NULL)
        return 0;
    return 1 + MSSample_length(self->next);
}

/// Sort by increasing age. Within age categories, nodes that are
/// inserted later sort as though they were older.
MSSample *MSSample_insert(MSSample *node, MSSample *new) {
    if(node == NULL)
        return new;
    int cmp = MSSample_compare(new, node);
    if(cmp >= 0) {
        node->next = MSSample_insert(node->next, new);
        return node;
    }
    new->next = node;
    return new;
}

/// Collapse adjacent nodes with the same age and population.
MSSample *MSSample_collapse(MSSample *node) {
    if(node == NULL)
        return NULL;
    while(node->next
          && node->age == node->next->age
          && node->pop == node->next->pop) {
        node->nsamp += node->next->nsamp;
        MSSample *tmp = node->next;
        node->next = node->next->next;
        free(tmp);
    }
    node->next = MSSample_collapse(node->next);
    return node;
}

/// Parse output of ms
void usage(void) {
    fprintf(stderr, "usage: ms2sim [options] [input_file_name]\n");
    fprintf(stderr, "   where options may include:\n");
    tellopt("-R <x> or --recombination <x>", "rate for adjacent nucleotides");
    tellopt("-h     or --help", "print this message");

    fprintf(stderr, " If no input file is given, ms2sim reads stdin.");
    fprintf(stderr, " ms2sim always writes\n  to standard output.\n");

    exit(EXIT_FAILURE);
}

// On input, glt should point to a tokenized string representing the
// ms command line, and nsamples should point to an int.
//
// The function returns a newly-allocated array of ints, whose dimension
// is *sampleDim, the number of populations specified on the ms command line.
// The i'th entry in this array is the haploid sample size of population i.
//
// On error, the function returns NULL.
unsigned *countSamples(GetLineTok *glt, int *sampleDim) {
    /*
      Here is the summary, from Hudson's msdoc.pdf, of the order in
      which haplotypes are printed in ms output: "The ancient DNA
      haplotypes are output after the haplotypes sampled in the
      present. nsam parameter is the number of chromosomes sampled at
      the present time. If multiple -eA switches are employed the
      ancient DNA haplotypes are output after the present day
      haplotypes in order from most recent to most ancient. If ancient
      DNA is sampled at the same time from two or more different
      subpopulations, the haplotypes are output in order of the -eA
      switches on the command line."

      In other words, samples are sorted first by age (youngest first)
      and then by the order in which they appear on the command
      line. The modern sample sizes are specified as arguments to -I
      and the ancient ones as arguments to -eA.
    */

    if(strcmp("ms", GetLineTok_token(glt, 0)) != 0) {
        fprintf(stderr,"%s:%d: input file is not ms output\n",
                __FILE__,__LINE__);
        return NULL;
    }

    int i, j, k, npops=0, h;
    char *token, *end;
    MSSample *mss = NULL;
    int ntokens = GetLineTok_ntokens(glt);

    // Read through tokens, looking for -I and -eA. Use these arguments
    // to set npops and push sample sizes into mss.
    i = 1;
    while(i < ntokens) {
        token = GetLineTok_token(glt, i);
        if(strcmp("-I", token) == 0) {
            assert(npops == 0);
            assert(mss == NULL);
            token = GetLineTok_token(glt, i+1);
            npops = strtol(token, &end, 10);
            if(end==token || npops <= 0) {
                fprintf(stderr,"%s:%d: illegal number of populations: %s\n",
                        __FILE__,__LINE__, token);
                goto fail;
            }
            for(j=0; j<npops; ++j) {
                token = GetLineTok_token(glt, i+j+2);
                h = strtol(token, &end, 10);
                if(end==token || h < 0) {
                    fprintf(stderr,"%s:%d: illegal number of populations: %s\n",
                            __FILE__,__LINE__, token);
                    goto fail;
                }
                // time=0.0, pop=j+1, sample size = h
                MSSample *new = MSSample_new(0.0, j+1, h);
                mss = MSSample_insert(mss, new);
            }
            i += npops + 2;
            continue;
        }else if(strcmp("-eA", token) == 0) {
            token = GetLineTok_token(glt, i+1);
            double t = strtod(token, &end);
            if(end==token || t < 0.0) {
                fprintf(stderr,"%s:%d: illegal time: %s\n",
                        __FILE__,__LINE__, token);
                goto fail;
            }
            token = GetLineTok_token(glt, i+2);
            int pop = strtol(token, &end, 10);
            if(end==token || pop < 1) {
                fprintf(stderr,"%s:%d: illegal population: %s\n",
                        __FILE__,__LINE__, token);
                goto fail;
            }
            token = GetLineTok_token(glt, i+3);
            h = strtol(token, &end, 10);
            if(end==token || h <= 0) {
                fprintf(stderr,"%s:%d: illegal sample size: %s\n",
                        __FILE__,__LINE__, token);
                goto fail;
            }
            // time=t, pop=pop, sample size = h
            MSSample *new = MSSample_new(t, pop, h);
            mss = MSSample_insert(mss, new);
            i += 4;
            continue;
        }
        i += 1;
    }

    mss = MSSample_collapse(mss);
    npops = MSSample_length(mss);

    *sampleDim = npops;
    unsigned *nsamples = malloc(npops * sizeof(nsamples[0]));
    CHECKMEM(nsamples);
    MSSample *node = mss;

    for(i=0; i < npops; ++i) {
        assert(node);
        nsamples[i] = node->nsamp;
        node = node->next;
    }
    MSSample_free(mss);

    return nsamples;

 fail:
    MSSample_free(mss);
    return NULL;
}

int main(int argc, char **argv) {
    long        maxtokens = 1<<8;
    long        inBuffSize = 1<<9;
    long        i, seq, site, ntokens, status;
    double      recombRate = 0.0;   /* rate for adjacent nucleotides */
    long        nseg, nseq;
    long       *nucpos;         /* vector of nucleotide positions */
    double     *mappos;         /* vector of map positions (cM)   */
    char       *buff=NULL, *buff2=NULL;
    int         optndx;
    time_t      currtime = time(NULL);

    /*
     * A matrix of dimension nseq X nseg. To access element for
     * given seq and site: m[nseg*seq + site];
     */
    char       *m;

    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
        {"recombination", required_argument, 0, 'R'},
        {"help", no_argument, 0, 'h'},
        {NULL, 0, NULL, 0}
    };

    printf("# ms2sim was run %s", ctime(&currtime));
    printf("# %-12s =", "cmd line");
    for(i = 0; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');

    /* flag arguments */
    for(;;) {
        i = getopt_long(argc, argv, "R:h", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            usage();
            break;
        case 'R':
            recombRate = strtod(optarg, NULL);
            break;
        case 'h':
            usage();
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

    // allocate memory
    GetLineTok *glt = GetLineTok(inBuffSize, maxtokens, stdin);
    CHECKMEM(glt);

    status = GetLineTok_next(glt, " \t", " \t\n");
    switch(status) {
    case EOF:
        DIE("No input");
        break;
    case ENOMEM:
        DIE("out of memory");
        break;
    default:
        DIE("this shouldn't happen");
        break;
    }

    // First line of input has MS command line
    int sampleDim=0;
    unsigned *nsamples = countSamples(glt, &sampleDim);

    if(n != sampleDim) {
        fprintf(stderr, 
                "%s:%d: number of labels on command line doesn't match\n"
                " number of subdivisions in ms input.\n",
                __FILE__,__LINE__);
        fprintf(stderr,"  %d labels != %d subdivisions\n",
                n, sampleDim);
        exit(EXIT_FAILURE);
    }

    printf("# %-12s =", "ms cmd");
    for(i = 0; i < ntokens; ++i)
        printf(" %s", GetLineTok_token(glt, (int) i));
    putchar('\n');
    fflush(stdout);

    nseq = strtol(GetLineTok_token(glt, 1), NULL, 10);
    nseg = -1;

    // read input through nseg assignment
    while(1) {
        status = GetLineTok_next(glt, " \t\n", " \t\n");
        switch(status) {
        case EOF:
            DIE("segsites not found in input");
            break;
        case ENOMEM:
            DIE("out of memory");
            break;
        default:
            DIE("this shouldn't happen");
            break;
        }

        ntokens = GetLineTok_ntokens(glt);
        if(ntokens == 0)
            continue;

        if(strcmp(Tokenizer_token(tkz, 0), "segsites") == 0) {
            // got segsites: assign and break out of loop
            nseg = strtol(Tokenizer_token(tkz, 1), NULL, 10);
            break;
        }
    }
    if(nseg <= 0) {
        fprintf(stderr, "%s:%d: Number of segregating sites is %ld",
                __FILE__, __LINE__, nseg);
        exit(EXIT_FAILURE);
    }

    status = GetLineTok_next(glt, ":", " \t\n");
    switch(status) {
    case EOF:
        DIE("positions not found in input");
        break;
    case ENOMEM:
        DIE("out of memory");
        break;
    default:
        DIE("this shouldn't happen");
        break;
    }

    ntokens = GetLineTok_ntokens(glt);
    if(ntokens != 2)
        DIE("positions not found in input");

    printf("# %-12s = %ld\n", "segsites", nseg);
    printf("# %-12s = %ld\n", "nsequences", nseq);
    printf("# %-12s = %d\n", "ploidy", 1);

    m = malloc(nseq * nseg * sizeof(m[0]));
    CHECKMEM(m);

    // Read data, putting values into array "m"
    buff = malloc((100 + nseg) * sizeof(buff[0]));
    CHECKMEM(buff);
    seq = 0;
    while(NULL != fgets(buff, inBuffSize, stdin)) {

        if(!strchr(buff, '\n') && !feof(stdin)) {
            fprintf(stderr, "%s:%d: input buffer overflow."
                    " Curr buff size: inBuffSize=%d\n",
                    __FILE__, __LINE__, inBuffSize);
            exit(EXIT_FAILURE);
        }

        // make sure input line has right number of sites
        i = strlen(buff);
        if(i == 0)
            continue;

        if(seq >= nseq) {
            fprintf(stderr,"%s:%d: seq=%ld >= nseq=%ld",
                   __FILE__, __LINE__, seq, nseq);
            exit(EXIT_FAILURE);
        }

        while(isspace(buff[i - 1])) {
            --i;
            buff[i] = '\0';
        }
        if(i != nseg) {
            fprintf(stderr,
                    "%s:%d input line has %ld chars; should have %d\n",
                    __FILE__, __LINE__, i, nseg);
            exit(EXIT_FAILURE);
        }
        memcpy(m + nseg * seq, buff, nseg * sizeof(m[0]));
        seq += 1;
    }
    assert(seq == nseq);

    // Output loop
    printf("#%9s %10s %14s %7s %s\n", "snp_id", "nucpos", "mappos",
           "alleles", "genotypes");
    for(site = 0; site < nseg; ++site) {
        printf("%10ld %10ld %14.12lf %7s ",
               site, nucpos[site], mappos[site], "01");

        for(seq = 0; seq < nseq; ++seq)
            printf("%c", m[nseg * seq + site]);
        putchar('\n');
    }

    GetLineTok_free(glt);
    free(nucpos);
    free(mappos);
    free(buff);
    free(buff2);
    free(m);

    return 0;
}
