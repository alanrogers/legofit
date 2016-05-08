/**
 * @file lego.c
 * @brief Simulate branch lengths
 *
 * @copyright Copyright (c) 2015, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "binary.h"
#include "branchtab.h"
#include "gptree.h"
#include "hashtab.h"
#include "jobqueue.h"
#include "misc.h"
#include "parse.h"
#include "patprob.h"
#include "parstore.h"
#include "sampndx.h"
#include <assert.h>
#include <float.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

char       *patLbl(size_t n, char buff[n], tipId_t tid, SampNdx * sndx);
void        usage(void);
int         comparePtrs(const void *void_x, const void *void_y);

void usage(void) {
    fprintf(stderr, "usage: lego [options] input_file_name\n");
    fprintf(stderr, "   where options may include:\n");
    tellopt("-i <x> or --nItr <x>", "number of iterations in simulation");
    tellopt("-t <x> or --threads <x>", "number of threads (default is auto)");
    tellopt("-h or --help", "print this message");
    exit(1);
}

/// Compare pointers to pointers to two tipId_t values.
///
/// @param void_x,void_y pointers to pointers to tipId_t values
/// @returns <0, 0, or >0 depending on whether the first arg is <,
/// ==, or > the second.
int comparePtrs(const void *void_x, const void *void_y) {
    tipId_t * const * x = (tipId_t * const *) void_x;
    tipId_t * const * y = (tipId_t * const *) void_y;

    // Major sort is on the number of samples
    // represented in the site pattern. Patterns with
    // fewer samples come first.
    int diff1bits = num1bits(**x) - num1bits(**y);
    if(diff1bits)
        return diff1bits;

    // Reverse order of bits so that low-order bit
    // is most significant. This ensures that the
    // sort order of samples corresponds to the
    // order in which they were listed in the input
    // data.
    unsigned rx = reverseBits(**x);
    unsigned ry = reverseBits(**y);

    return ry - rx;
}

/// Generate a label for site pattern tid. Label goes into
/// buff. Function returns a pointer to buff;
char       *patLbl(size_t n, char buff[n], tipId_t tid, SampNdx * sndx) {
    int         maxbits = 40;
    int         bit[maxbits];
    int         i, nbits;
    nbits = getBits(tid, maxbits, bit);
    buff[0] = '\0';
    char        lbl[100];
    for(i = 0; i < nbits; ++i) {
        snprintf(lbl, sizeof(lbl), "%s",
                 SampNdx_lbl(sndx, (unsigned) bit[i]));
        if(strlen(buff) + strlen(lbl) >= n)
            eprintf("%s:%s:%d: buffer overflow\n", __FILE__, __func__,
                    __LINE__);
        strcat(buff, lbl);
        if(i + 1 < nbits && 1 + strlen(buff) < n)
            strcat(buff, ":");
    }
    return buff;
}

int main(int argc, char **argv) {

    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
        {"nItr", required_argument, 0, 'i'},
        {"threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {NULL, 0, NULL, 0}
    };

    printf("#################################################\n"
           "# lego: estimate probabilities of site patterns #\n"
           "#################################################\n");
    putchar('\n');

    int         i, j;
    time_t      currtime = time(NULL);
	double      lo_twoN = 0.0, hi_twoN = 1e6;  // twoN bounds
	double      lo_t = 0.0, hi_t = 1e6;        // t bounds
#if defined(__DATE__) && defined(__TIME__)
    printf("# Program was compiled: %s %s\n", __DATE__, __TIME__);
#endif
    printf("# Program was run: %s\n", ctime(&currtime));

    printf("# cmd:");
    for(i = 0; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');

    fflush(stdout);

    int         nTasks = 0;     // total number of tasks
    int         optndx;
    unsigned long nreps = 100;
    char        fname[200] = { '\0' };

    // command line arguments
    for(;;) {
        i = getopt_long(argc, argv, "i:t:h", myopts, &optndx);
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
        case 't':
            nTasks = strtol(optarg, NULL, 10);
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

    if(nTasks == 0)
        nTasks = getNumCores();

    if(nTasks > nreps)
        nTasks = nreps;

    // Divide repetitions among tasks.
    long        reps[nTasks];
    {
        ldiv_t      qr = ldiv((long) nreps, (long) nTasks);
        assert(qr.quot > 0);
        for(j = 0; j < nTasks; ++j)
            reps[j] = qr.quot;
        assert(qr.rem < nTasks);
        for(j=0; j < qr.rem; ++j)
            reps[j] += 1;
#ifndef NDEBUG
        // make sure the total number of repetitions is nreps.
        long        sumreps = 0;
        for(j = 0; j < nTasks; ++j) {
            assert(reps[j] > 0);
            sumreps += reps[j];
        }
        assert(sumreps = nreps);
#endif
    }

	int maxpat = 10;
	tipId_t pat[maxpat];
	double prob[maxpat];
	Bounds bnd = {
			.lo_twoN = lo_twoN,
			.hi_twoN = hi_twoN,
			.lo_t = lo_t,
			.hi_t = hi_t
	};

    printf("# nreps       : %lu\n", nreps);
    printf("# nthreads    : %d\n", nTasks);
    printf("# input file  : %s\n", fname);

	unsigned npat = patprob(maxpat, pat, prob, nTasks, reps, 0, fname, bnd);

    // Determine order for printing lines of output
    tipId_t  *ptr[npat];
    unsigned ord[npat];
    for(j=0; j < npat; ++j)
        ptr[j] = pat+j;
    qsort(ptr, (size_t) npat, sizeof(ptr[0]), comparePtrs);
    for(i=0; i<npat; ++i)
        ord[i] = ptr[i]-pat;

    printf("#%14s %10s\n", "SitePat", "Prob");
    char        buff[100];
    for(j = 0; j < npat; ++j) {
        char        buff2[100];
        snprintf(buff2, sizeof(buff2), "%s",
                 patLbl(sizeof(buff), buff, pat[ord[j]], &(taskarg[0]->sndx)));
        printf("%15s %10.7lf\n", buff2, branchLength[ord[j]]);
    }

    return 0;
}
