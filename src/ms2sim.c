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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <getopt.h>
#include <assert.h>
#include <math.h>
#include "misc.h"
#include "ini.h"
#include "tokenizer.h"

void        usage(void);

/**
 * Parse output of MS, produce gtp format, the input for eld.
 */
void usage(void) {
    fprintf(stderr, "usage: ms2sim [options] [input_file_name]\n");
    fprintf(stderr, "   where options may include:\n");
    tellopt("-t <x> or --maxTokens <x>", "maximum tokens per input line");
    tellopt("-b <x> or --inBuffSize <x>", "size of input buffer");
    tellopt("-R <x> or --recombination <x>", "rate for adjacent nucleotides");
    tellopt("-h     or --help", "print this message");

    fprintf(stderr, "  To figure out appropriate values for -t and -b,"
            " execute the following at\n  the command line:\n\n"
            "     grep positions ms.out | wc\n\n"
            "  where \"ms.out\" is the ms output file.");
    fprintf(stderr, " This command prints 3 numbers.\n"
            "  Number of tokens must be as large as 2nd;"
            " buffer size as large as the\n  3rd.");
    fprintf(stderr, " If no input file is given, ms2sim reads stdin.");
    fprintf(stderr, " ms2sim always writes\n  to standard output.\n");

    exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
    long        maxtokens = 1000;
    long        inBuffSize = 1000;
    long        i, seq, site, ntokens;
    double      recombRate = 0.0;   /* rate for adjacent nucleotides */
    double      twoN0 = 0.0;    /* pop size in epoch 0 */
    char       *ifname = NULL;
    FILE       *ifp = NULL;
    Tokenizer  *tkz = NULL;
    long        nseg, nseq, nnuc;
    long       *nucpos;         /* vector of nucleotide positions */
    double     *mappos;         /* vector of map positions (cM)   */
    char       *buff, *buff2;
    int         optndx;
    time_t      currtime = time(NULL);

    /*
     * A matrix of dimension nseq X nseg. To access element for
     * given seq and site: m[nseg*seq + site];
     */
    char       *m;

    /* import definitions from initialization file */
    Ini        *ini = Ini_new(INIFILE);

    if(ini) {
        Ini_setDbl(ini, "recombination", &recombRate, !MANDATORY);
        twoN0 = Ini_twoN0(ini);
        Ini_free(ini);
        ini = NULL;
    }

    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
        {"maxTokens", required_argument, 0, 't'},
        {"inBuffSize", required_argument, 0, 'b'},
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
        i = getopt_long(argc, argv, "t:b:", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            usage();
            break;
        case 't':
            maxtokens = strtol(optarg, NULL, 10);
            break;
        case 'b':
            inBuffSize = strtol(optarg, NULL, 10);
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
        if(ifp == NULL)
            eprintf("ERR@%s:%d: couldn't open file \"%s\"\n",
                    __FILE__, __LINE__, ifname);
        break;
    default:
        fprintf(stderr, "Only one input file is allowed\n");
        usage();
    }
    assert(ifname != NULL);

    if(ifp == NULL)
        usage();

    /* allocate memory */
    tkz = Tokenizer_new(maxtokens);
    buff = malloc(inBuffSize * sizeof(buff[0]));
    checkmem(buff, __FILE__, __LINE__);

    do {
        if(NULL == fgets(buff, inBuffSize, ifp))
            eprintf("ERR@%s:%d: Input file \"%s\" is empty\n",
                    __FILE__, __LINE__, ifname);
    } while(strempty(buff) || strcomment(buff));

    /* First line of input has MS command line */
    if(!strchr(buff, '\n') && !feof(ifp))
        eprintf("ERR@%s:%d: input buffer overflow."
                " Curr buff size: inBuffSize=%d\n",
                __FILE__, __LINE__, inBuffSize);
    Tokenizer_split(tkz, buff, " \t\n");
    ntokens = Tokenizer_strip(tkz, " \t\n");

    if(ntokens < 1 || strcmp(Tokenizer_token(tkz, 0), "ms") != 0) {
        fprintf(stderr, "1st line looks wrong:\n ");
        int         truncated = 0;

        if(ntokens > 10) {
            truncated = ntokens - 10;
            ntokens = 10;
        }
        for(i = 0; i < ntokens; ++i)
            fprintf(stderr, " %s", Tokenizer_token(tkz, i));
        if(truncated)
            fprintf(stderr, " <%d more tokens>", truncated);
        putchar('\n');
        exit(1);
    }

    printf("# %-12s =", "sim cmd");
    for(i = 0; i < ntokens; ++i)
        printf(" %s", Tokenizer_token(tkz, (int) i));
    putchar('\n');
    fflush(stdout);

    i = Tokenizer_find(tkz, "-r");
    if(i < ntokens) {
        double      fourNcMax = 0.0;

        ++i;
        if(i >= ntokens)
            eprintf("ERR@%s:%d: missing argument to -r in ms command line",
                    __FILE__, __LINE__);
        fourNcMax = strtod(Tokenizer_token(tkz, (int) i), NULL);

        ++i;
        if(i >= ntokens)
            eprintf("ERR@%s:%d: missing argument to -r in ms command line",
                    __FILE__, __LINE__);
        nnuc = strtol(Tokenizer_token(tkz, (int) i), NULL, 10);

        /* Make sure recombination rate makes sense */
        if(twoN0 > 0.0) {
            double      fourNcMaxExpected =
                recombRate * twoN0 * 2 * (nnuc - 1);
            double      err = fabs(fourNcMax - fourNcMaxExpected) / fourNcMax;

            if(err > 1e-4) {
                eprintf("ERR@%s:%d: fourNcMax=%lg; should be %lg; relerr=%le",
                        __FILE__, __LINE__,
                        fourNcMax, fourNcMaxExpected, err);
                exit(1);
            }
        }
    } else {
        fprintf(stderr, "ERR: -r absent from ms command line\n");
        exit(EXIT_FAILURE);
    }

    nseq = strtol(Tokenizer_token(tkz, 1), NULL, 10);
    nseg = -1;

    /* read input through nseg assignment */
    while(1) {

        if(NULL == fgets(buff, inBuffSize, ifp)) {
            fprintf(stderr, "Error: segsites not found in input\n");
            exit(1);
        }

        if(!strchr(buff, '\n') && !feof(ifp))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, inBuffSize);

        Tokenizer_split(tkz, buff, ":");
        ntokens = Tokenizer_strip(tkz, " \t\n");
        if(ntokens == 0)
            continue;

        if(strcmp(Tokenizer_token(tkz, 0), "segsites") == 0) {
            /* got segsites: assign and break out of loop */
            nseg = strtol(Tokenizer_token(tkz, 1), NULL, 10);
            break;
        }
    }
    if(nseg <= 0)
        eprintf("ERR@%s:%d: Number of segregating sites is %ld",
                __FILE__, __LINE__, nseg);
    if(nseg > maxtokens)
        eprintf("ERR@%s:%d: The data have %ld segregating sites,"
                " but I can accomodate only %ld.\n"
                "Rerun using \"--maxTokens %ld\".\n",
                __FILE__, __LINE__, nseg, maxtokens, nseg + 1);

    if(NULL == fgets(buff, inBuffSize, ifp)) {
        fprintf(stderr, "Error: positions not found in input\n");
        exit(1);
    }

    if(!strchr(buff, '\n') && !feof(ifp))
        eprintf("ERR@%s:%d: input buffer overflow."
                " buff size: %d", __FILE__, __LINE__, inBuffSize);

    Tokenizer_split(tkz, buff, ":");
    ntokens = Tokenizer_strip(tkz, " \t\n");
    if(ntokens != 2)
        eprintf("ERR@%s:%d: positions not found in input",
                __FILE__, __LINE__);

    buff2 = strdup(Tokenizer_token(tkz, 1));
    Tokenizer_split(tkz, buff2, " \t");
    ntokens = Tokenizer_strip(tkz, " \t");
    if(ntokens != nseg) {
        fflush(stdout);
        fprintf(stderr, "@%s:%d:Parsing positions: ntokens=%ld != nseg=%ld\n",
                __FILE__, __LINE__, ntokens, nseg);
        Tokenizer_printSummary(tkz, stderr);
        exit(EXIT_FAILURE);
    }

    nucpos = malloc(nseg * sizeof(nucpos[0]));
    checkmem(nucpos, __FILE__, __LINE__);

    mappos = malloc(nseg * sizeof(mappos[0]));
    checkmem(mappos, __FILE__, __LINE__);

    for(i = 0; i < ntokens; ++i) {
        double      pos = strtod(Tokenizer_token(tkz, (int) i), 0);

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
    checkmem(m, __FILE__, __LINE__);

    /* Read data, putting values into array "m" */
    seq = 0;
    while(NULL != fgets(buff, inBuffSize, ifp)) {

        if(!strchr(buff, '\n') && !feof(ifp))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " Curr buff size: inBuffSize=%d\n",
                    __FILE__, __LINE__, inBuffSize);

        /* make sure input line has right number of sites */
        i = strlen(buff);
        if(i == 0)
            continue;

        if(seq >= nseq) {
            printf("ERR@%s:%d: seq=%ld >= nseq=%ld",
                   __FILE__, __LINE__, seq, nseq);
            exit(EXIT_FAILURE);
        }

        while(isspace(buff[i - 1])) {
            --i;
            buff[i] = '\0';
        }
        if(i != nseg)
            eprintf("ERR@%s:%d input line has %ld chars; should have %d\n",
                    __FILE__, __LINE__, i, nseg);
        memcpy(m + nseg * seq, buff, nseg * sizeof(m[0]));
        seq += 1;
    }
    assert(seq == nseq);

    /* Output loop */
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
    if(ifname)
        free(ifname);
    free(nucpos);
    free(mappos);
    free(buff);
    free(m);
    free(buff2);

    return 0;
}
