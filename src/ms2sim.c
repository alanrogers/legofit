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
 *     # nnucleotides = 1000000
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
    int nsamp;    // number of haploid samples
    MSSample *next;
};

void      usage(void);
unsigned *countSamples(Tokenizer *tkz, int *sampleDim);
MSSample *MSSample_new(double age, int nsamp);
int       MSSample_compare(MSSample *left, MSSample *right);
MSSample *MSSample_insert(MSSample *node, MSSample *new);

MSSample *MSSample_new(double age, int position, int nsamp) {
    MSSample *self = malloc(sizeof(MSSample));
    CHECKMEM(self);
    memset(self, 0, sizeof(MSSample));
    return self;
}

// Sort first by age then by position.
int MSSample_compare(MSSample *left, MSSample *right) {
    if(left->age > right->age)
        return 1;
    if(left->age < right->age)
        return -1;
    return 0;
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

/**
 * Parse output of MS, produce gtp format, the input for eld.
 */
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

    int i, j, k;
    long h;
    char *token, *end;
    UINTqueue **queue=NULL;  // array of UINTqueue objects, one per population
    int ntokens = GetLineTok_ntokens(glt);
    int npops=0;

    // Read through tokens, looking for -I and -eA. Use these arguments
    // to set npops and array of queues.
    for(i=1; i < ntokens; ++i) {
        token = GetLineTok_token(glt, i);
        if(strcmp("-I", token) == 0) {
            if(npops == 0) {
                // count populations and allocate nsamples
                assert(queue == NULL);
                for(j=i+2; j<ntokens; ++j) {
                    token = GetLineTok_token(glt, j);
                    h = strtol(token, &end, 10);
                    if(end==token || h < 0) // token isn't a nonnegative int
                        break;
                    else           // token is a nonnegative int
                        ++npops;
                }
                if(npops == 0) {
                    fprintf(stderr,"%s:%d: npops is zero\n",
                            __FILE__,__LINE__);
                    goto fail:
                }
                queue = malloc(npops * sizeof(queue[0]));
                CHECKMEM(queue);
                memset(queue, 0, npops * sizeof(queue[0]));
            }
            // increment queues
            assert(npops > 0);
            assert(queue != NULL);
            for(j=0; j < npops; ++j) {
                token = GetLineTok_token(glt, i+2+j);
                h = strtol(token, &end, 10);
                if(end == token || h < 0) {
                    fprintf(stderr,"%s:%d: read \"%s\" when"
                            " expecting a sample size\n",
                            __FILE__, __LINE__, token);
                    goto fail;
                }
                if(h>0)
                    queue[j] = UINTqueue_push(queue[j], (unsigned) h);
            }
            // advance to last argument of -I or -eI
            i += 1 + npops;
        }else if(strcmp("-eA", token) == 0) {
            token = GetLineTok_token(glt, i+2); // index of pop
            j = -1 + strtod(token, &end, 10);
            if(end == token || j < 0) {
                fprintf(stderr,"%s:%d: illegal population index %d\n",
                        __FILE__,__LINE__, j);
                goto fail;
            }
            if(j >= npops) {
                fprintf(stderr,"%s:%d: population index (%d) >="
                        " npops (%d)\n",
                        __FILE__,__LINE__, j, npops);
                goto fail;
            }
            token = GetLineTok_token(glt, i+3); // index of pop
            h = strtod(token, &end, 10);
            if(end == token || h < 0) {
                fprintf(stderr,"%s:%d: read \"%s\" when"
                        " expecting a sample size\n",
                        __FILE__, __LINE__, token);
                goto fail;
            }
        }
    }

    *sampleDim=0;
    for(j=0; j < npops; ++j)
        *sampleDim += UINTqueue_length(queue[j]);
    assert(*sampleDim > 0);
    unsigned *nsamples = malloc(*sampleDim * sizeof(nsamples[0]));
    CHECKMEM(nsamples);
    for(i=j=0; i < npops; ++i) {
        unsigned n;
        while(queue[i]) {
            queue[i] = UINTqueue_pop(queue[i], &n);
            nsamples[j++] = n;
        }
    }
    assert(j == *sampleDim);
#ifndef NDEBUG
    for(i=0; i < npops; ++i)
        assert(queue[i] == NULL);
#endif
    free(queue);

    return nsamples;

 fail:
    if(queue != NULL) {
        for(k=0; k<npops; ++k) {
            if(queue[k] != NULL)
                queue[k] = UINTqueue_free(queue[k]);
        }
        free(queue);
    }
    return NULL;
}

int main(int argc, char **argv) {
    long        maxtokens = 1<<8;
    long        inBuffSize = 1<<9;
    long        i, seq, site, ntokens, status;
    double      recombRate = 0.0;   /* rate for adjacent nucleotides */
    char       *ifname = NULL;
    FILE       *ifp = NULL;
    Tokenizer  *tkz = NULL;
    long        nseg, nseq, nnuc;
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

    /* remaining option gives file name */
    switch (argc - optind) {
    case 0:
        ifname = strdup("stdin");
        ifp = stdin;
        fprintf(stderr, "Reading from standard input\n");
        break;
    case 1:
        ifname = strdup(argv[optind]);
        ifp = fopen(ifname, "r");
        if(ifp == NULL) {
            fprintf(stderr, "%s:%d: couldn't open file \"%s\"\n",
                    __FILE__, __LINE__, ifname);
            exit(EXIT_FAILURE);
        }
        break;
    default:
        fprintf(stderr, "Only one input file is allowed\n");
        usage();
    }
    assert(ifname != NULL);

    if(ifp == NULL)
        usage();

    // allocate memory
    GetLineTok *glt = GetLineTok(inBuffSize, maxtokens, ifp);
    CHECKMEM(glt);

    status = GetLineTok_next(glt, " \t\n", " \t\n");
    switch(status) {
    case EOF:
        DIE("Can't read input file");
        break;
    case ENOMEM:
        DIE("out of memory");
        break;
    default:
        DIE("this shouldn't happen");
        break;
    }

    // First line of input has MS command line
    ntokens = GetLineTok_ntokens(glt);

    if(ntokens < 1 || strcmp(GetLineTok_token(glt, 0), "ms") != 0) {
        fprintf(stderr, "1st line looks wrong:\n");
        int         truncated = 0;

        if(ntokens > 10) {
            truncated = ntokens - 10;
            ntokens = 10;
        }
        for(i = 0; i < ntokens; ++i)
            fprintf(stderr, " %s", GetLineTok_token(glt, i));
        if(truncated)
            fprintf(stderr, " <%d more tokens>", truncated);
        putchar('\n');
        exit(EXIT_FAILURE);
    }

    printf("# %-12s =", "ms cmd");
    for(i = 0; i < ntokens; ++i)
        printf(" %s", GetLineTok_token(glt, (int) i));
    putchar('\n');
    fflush(stdout);

    i = GetLineTok_find(glt, "-r");
    if(i < ntokens) {
        double      fourNcMax = 0.0;

        ++i;
        if(i >= ntokens)
            DIE("missing argument to -r in ms command line");
        fourNcMax = strtod(GetLineTok_token(glt, (int) i), NULL);

        ++i;
        if(i >= ntokens)
            DIE("missing argument to -r in ms command line");
        nnuc = strtol(GetLineTok_token(glt, (int) i), NULL, 10);
    } else
        DIE("-r absent from ms command line");

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

    Tokenizer *tkz = Tokenizer_new(0);
    buff2 = strdup(GetLineTok_token(glt, 1));
    Tokenizer_split(tkz, buff2, " \t");
    ntokens = Tokenizer_strip(tkz, " \t");
    if(ntokens != nseg) {
        fflush(stdout);
        fprintf(stderr, "%s:%d:Parsing positions: ntokens=%ld != nseg=%ld\n",
                __FILE__, __LINE__, ntokens, nseg);
        Tokenizer_printSummary(tkz, stderr);
        exit(EXIT_FAILURE);
    }

    nucpos = malloc(nseg * sizeof(nucpos[0]));
    CHECKMEM(nucpos);

    mappos = malloc(nseg * sizeof(mappos[0]));
    CHECKMEM(mappos);

    for(i = 0; i < ntokens; ++i) {
        double pos = strtod(Tokenizer_token(tkz, (int) i), 0);

        pos *= nnuc;            /* proportion of sites --> num of sites */
        mappos[i] = pos * recombRate * 0.01;
        pos = floor(pos + 0.5); /* round to int */
        nucpos[i] = (long int) pos;
    }

    printf("# %-12s = %ld\n", "nnucleotides", nnuc);
    printf("# %-12s = %ld\n", "segsites", nseg);
    printf("# %-12s = %ld\n", "nsequences", nseq);
    printf("# %-12s = %lg\n", "recomb rate per site", recombRate);
    printf("# %-12s = %d\n", "ploidy", 1);

    m = malloc(nseq * nseg * sizeof(m[0]));
    CHECKMEM(m);

    // Read data, putting values into array "m"
    buff = malloc((100 + nseg) * sizeof(buff[0]));
    CHECKMEM(buff);
    seq = 0;
    while(NULL != fgets(buff, inBuffSize, ifp)) {

        if(!strchr(buff, '\n') && !feof(ifp)) {
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

    fprintf(stderr, "Maximum mappos: %lf\n", mappos[site - 1]);
    fprintf(stderr, "Should be     : %lf\n", nnuc * recombRate * 0.01);

    fclose(ifp);
    Tokenizer_free(tkz);
    GetLineTok_free(glt);
    if(ifname)
        free(ifname);
    free(nucpos);
    free(mappos);
    free(buff);
    free(buff2);
    free(m);

    return 0;
}
