/**
@file raf2daf.c
@page raf2daf
@brief Convert raf files to daf files

# raf2daf: converts raf files to daf files

Raf2daf reads data in .raf format and writes the corresponding daf files. 

# Usage

    Usage: raf2daf [options] <in1> <in2> ...
       where <in1> and <in2> are input files in raf format. The last
       input file should be the outgroup. Writes to standard
       output.

    Options may include:
       -m or --logMismatch
          log AA/DA mismatches to raf2daf.log
       -A or --logAA
          log sites with uncallable ancestral alleles
       -a or --logAll
          log all sites to raf2daf.log
       -h or --help
          Print this message

@copyright Copyright (c) 2018, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "rafreader.h"
#include "misc.h"
#include "strint.h"
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

static void usage(void);

const char *useMsg =
    "\nUsage: raf2daf [options] <in_1> <in_2> <in_3>... \n"
    "   where <in_i> are input files in raf format, the last of which\n"
    "   should be the outgroup. Output daf files have names like those\n"
    "   of input files but with raf changed to daf. At least 3 input\n"
    "   files are required";

/// Print usage message and die.
static void usage(void) {
    fputs(useMsg, stderr);
    fputs("\nOptions may include:\n", stderr);
    tellopt("-m or --logMismatch", "Log REF mismatches to raf2daf.log");
    tellopt("-A or --logAA", "Log sites with uncallable ancestral allele");
    tellopt("--version", "Print version and exit");
    tellopt("-h or --help", "Print this message");
    exit(1);
}

int main(int argc, char **argv) {
    int         i, j, status, optndx, done;
    char        errbuff[1024] = { '\0' };
    const char *logfname = "raf2daf.log";
    int         logMismatch = 0, logAA = 0;
    FILE       *logfile = NULL;

    static struct option myopts[] = {
        // {char *name, int has_arg, int *flag, int val}
        {"logMismatch", no_argument, 0, 'm'},
        {"logAA", no_argument, 0, 'A'},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'V'},
        {NULL, 0, NULL, 0}
    };

    // command line arguments
    for(;;) {
        i = getopt_long(argc, argv, "hmAV", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case 'V':
            printf("raf2daf version %s\n", GIT_VERSION);
            return 0;
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
    if(m < 2) {
        fprintf(stderr,"At least 3 input files are required\n");
        usage();
    }

    char       *ifname[n];
    FILE       *ofp[m];
    RAFReader  *r[n];

    // Parse remaining arguments, each of which should be the name of
    // an input file.
    for(i = 0; i < n; ++i) {
        ifname[i] = argv[i + optind];
        r[i] = RAFReader_new(ifname[i]);
        if(i == n-1)
            continue;
        char *start = strrchr(ifname[i], '/');
        if(start == NULL)
            start = ifname[i];
        else
            ++start;
        char *ofname = strdup(start);
        int k = strlen(ofname);
        if(0 != strcmp("raf", ofname + k-3)) {
            fprintf(stderr,"%s:%d: input files should end with \"raf\".\n",
                    __FILE__,__LINE__);
            fprintf(stderr,"   got \"%s\".\n", ofname);
            exit(EXIT_FAILURE);
        }
        ofname[k-3] = 'd';
        ofp[i] = fopen(ofname, "w");
        if(ofp[i] == NULL) {
            fprintf(stderr,"%s:%d: can't open \"%s\" for output.\n",
                    __FILE__,__LINE__, ofname);
            exit(EXIT_FAILURE);
        }
        fprintf(stderr,"writing to %s\n", ofname);
        fprintf(ofp[i], "#%3s %10s %2s %2s %s\n",
                "chr", "pos", "aa", "da", "daf");
    }

    if(logMismatch || logAA) {
        logfile = fopen(logfname, "w");
        if(logfile == NULL) {
            fprintf(stderr, "Can't write to file \"%s\".\n", logfname);
            exit(EXIT_FAILURE);
        }
    }

    printf("raf2daf version %s\n", GIT_VERSION);

    unsigned long nsites = 0, nbadaa = 0, nbadref=0, nmultalt=0;

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

        for(j = 0; j < m; ++j) {
            // frequencies of reference and derived alleles
            double raf = RAFReader_raf(r[j]);
            double daf = RAFReader_daf(r[j]);
            const char *ancestral, *derived;
            if(raf == daf) {
                derived = RAFReader_ref(r[j]);
                ancestral = RAFReader_alt(r[j]);
            }else{
                derived = RAFReader_alt(r[j]);
                ancestral = RAFReader_ref(r[j]);
            }
            fprintf(ofp[j], "%4s %10lu %2s %2s %0.18g\n",
                    RAFReader_chr(r[j]),
                    RAFReader_nucpos(r[j]),
                    ancestral, derived, daf);
        }
    }
    printf("Aligned sites                  : %lu\n", nsites);
    if(nbadref)
        printf("Disagreements about ref allele : %lu\n", nbadref);
    if(nmultalt)
        printf("Sites with multiple alt alleles: %lu\n", nmultalt);
    if(nbadaa)
        printf("Undetermined ancestral allele  : %lu\n", nbadaa);
    printf("Sites used                     : %lu\n",
           nsites - nbadaa - nbadref - nmultalt);

    for(i = 0; i < n; ++i)
        RAFReader_free(r[i]);
    for(i=0; i < m; ++i)
        fclose(ofp[i]);
    if(logfile)
        fclose(logfile);
    fprintf(stderr, "raf2daf is finished\n");
    return 0;
}
