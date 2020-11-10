/**
@file rafdaf.c
@page rafdaf
@brief Read a list of .raf  files and write derived allele frequencies.

# Rafdaf: read .raf files and write derived allele frequencies

Rafdaf reads multiple .raf files, the last of which must be an
outgroup, which is used to call ancestral alleles. It writes several 
columns to standard output. The first output column is "chr", a
character string representing the chromosome; the second is "pos", an
integer representing the nucleotide position along the chromosome. The
remaining columns give derived allele frequencies of the files listed
on the command line. The outgroup (which must appear last on the
command line) is not included in the output.

# Usage

    Usage: rafdaf [options] <x>=<in_1> <y>=<in_2> ... outgroup=<in_K>
       where <x> and <y> are arbitrary labels, and <in_i> are input
       files in raf format. Labels may not include the character ":".
       Final label must be "outgroup". Writes to standard output.

       If input file name ends with .gz, input is decompressed using
       gunzip.

    Options may include:
       --version
          Print version and exit
       -h or --help
          Print this message

@copyright Copyright (c) 2019 Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "binary.h"
#include "rafreader.h"
#include "misc.h"
#include "strint.h"
#include "error.h"
#include "typedefs.h"
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void usage(void);

const char *useMsg =
    "\nUsage: rafdaf [options] <x>=<in_1> <y>=<in_2> ... outgroup=<in_K>\n"
    "   where <x> and <y> are arbitrary labels, and <in_i> are input\n"
    "   files in raf format. Labels may not include the character \":\".\n"
    "   Final label must be \"outgroup\". Writes derived allele frequencies\n"
    "   to standard output. Minimum number of input files: 2 plus outgroup.\n"
    "\n"
    "   If input file name ends with .gz, input is decompressed using\n"
    "   gunzip.\n";

/// Print usage message and die.
static void usage(void) {
    fputs(useMsg, stderr);
    fputs("\nOptions may include:\n", stderr);
    tellopt("--version", "Print version and exit");
    tellopt("-h or --help", "Print this message");
    exit(1);
}

int main(int argc, char **argv) {
    int         i, j, status, optndx, done;
    StrInt     *strint = StrInt_new();
    static struct option myopts[] = {
        // {char *name, int has_arg, int *flag, int val}
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'V'},
        {NULL, 0, NULL, 0}
    };
    char errbuff[100] = { '\0' };

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
        case 'V':
            printf("rafdaf version %s\n", GIT_VERSION);
            return 0;
        case 'h':
            usage();
            break;
        default:
            usage();
        }
    }

    // remaining options: input files
    int         n = argc - optind;  // number of input files
    int         m = n-1;            // number excluding outgroup
    if(m < 2)
        usage();

    char       *poplbl[n];
    char       *fname[n];
    LblNdx      lndx;
    LblNdx_init(&lndx);
    RAFReader  *r[n];

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

    printf("# rafdaf version %s\n", VERSION);
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

    unsigned long nsites = 0, nbadaa = 0, nbadref=0, nmultalt=0;

    // Header
    printf("chr pos");
    for(j=0; j<m; ++j)
        printf(" %s", poplbl[j]);
    putchar('\n');

    // Iterate through raf files
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
            continue;
        case MULTIPLE_ALT:
            ++nsites;
            ++nmultalt;
            continue;
        case NO_ANCESTRAL_ALLELE:
            ++nsites;
            ++nbadaa;
            continue;
        default:
            // something wrong.
            mystrerror_r(status, errbuff, sizeof errbuff);
            fprintf(stderr,"%s:%d: input error (%s)\n",
                    __FILE__,__LINE__, errbuff);
            exit(EXIT_FAILURE);
        }

        // Output
        printf("%s %lu", RAFReader_chr(r[0]), RAFReader_nucpos(r[0]));
        for(j = 0; j < m; ++j)
            printf(" %0.8lg", RAFReader_daf(r[j]));
        putchar('\n');
    }
    fprintf(stderr, "# Aligned sites                  : %lu\n", nsites);
    if(nbadref)
        fprintf(stderr,"# Disagreements about ref allele : %lu\n", nbadref);
    if(nmultalt)
        fprintf(stderr,"# Sites with multiple alt alleles: %lu\n", nmultalt);
    if(nbadaa)
        fprintf(stderr,"# Undetermined ancestral allele  : %lu\n", nbadaa);
    fprintf(stderr,"# Sites used                     : %lu\n",
           nsites - nbadaa - nbadref - nmultalt);

    for(i = 0; i < n; ++i)
        RAFReader_free(r[i]);
    StrInt_free(strint);
    return 0;
}
