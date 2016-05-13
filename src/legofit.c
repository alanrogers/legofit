/**
 * @file legofit.c
 * @brief Estimate parameters describing population sizes, the times of separations
 * and of episodes of gene flow, and levels of gene flow.
 *
 * @copyright Copyright (c) 2016, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "gptree.h"
#include "patprob.h"
#include "parstore.h"
#include "lblndx.h"
#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <sys/types.h>
#include <time.h>

void        usage(void);

void usage(void) {
    fprintf(stderr, "usage: legofit [options] input_file_name\n");
    fprintf(stderr, "   where options may include:\n");
    tellopt("-i <x> or --nItr <x>", "number of iterations in simulation");
    tellopt("-t <x> or --threads <x>", "number of threads (default is auto)");
    tellopt("-h or --help", "print this message");
    exit(1);
}


int main(int argc, char **argv) {

    static struct option myopts[] = {
        /* {char *name, int has_arg, int *flag, int val} */
        {"nItr", required_argument, 0, 'i'},
        {"threads", required_argument, 0, 't'},
		{"crossover", required_argument, 0, 'x'},
		{"scaleFactor", required_argument, 0, 'F'},
        {"strategy", required_argument, 0, 's'},
        {"deTol", required_argument, 0, 'V'},
        {"help", no_argument, 0, 'h'},
        {NULL, 0, NULL, 0}
    };

    printf("########################################\n"
           "# legofit: estimate population history #\n"
           "########################################\n");
    putchar('\n');

    int         i, j;
    time_t      currtime = time(NULL);
    double      lo_twoN = 0.0, hi_twoN = 1e6;  // twoN bounds
    double      lo_t = 0.0, hi_t = 1e6;        // t bounds

	// DiffEv parameters
	double      F = 0.9;
	double      CR = 0.8;
	double      deTol = 1e-4;
	int         nPts;
	int         strategy = 1;
	int         ptsPerDim = 10;

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
        i = getopt_long(argc, argv, "i:t:F:p:s:V:x:h", myopts, &optndx);
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
		case 'F':
			F = strtod(optarg, 0);
			break;
        case 'p':
            ptsPerDim = strtol(optarg, NULL, 10);
            break;
		case 's':
			strategy = strtol(optarg, NULL, 10);
            break;
        case 'V':
            deTol = strtod(optarg, 0);
            break;
		case 'x':
			CR = strtod(optarg, 0);
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

    // parameters for Differential Evolution
    DiffEvPar   dep = {
        .dim = phdim,
        .ptsPerDim = ptsPerDim,
        .genmax = nItr,
        .refresh = 3,  // how often to print a line of output
        .strategy = strategy,
        .nthreads = nthreads,
        .verbose = verbose,
        .seed = (unsigned long) time(NULL),
        .F = F,
        .CR = CR,
        .deTol = deTol,
        .loBound = loBnd,
        .hiBound = hiBnd,
		.jobData = &costPar,
        .objfun = costFun,
		.threadData = &tdata,
		.ThreadState_new = ThreadState_new,
		.ThreadState_free = ThreadState_free
    };

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

    printf("# nreps       : %lu\n", nreps);
    printf("# nthreads    : %d\n", nTasks);
    printf("# input file  : %s\n", fname);

    int maxpat = 10;
    tipId_t pat[maxpat];
    double prob[maxpat];
    Bounds bnd = {
            .lo_twoN = lo_twoN,
            .hi_twoN = hi_twoN,
            .lo_t = lo_t,
            .hi_t = hi_t
    };
    GPTree *gptree = GPTree_new(fname, bnd);
	LblNdx lblndx;

    unsigned npat = patprob(maxpat, pat, prob, gptree, &lblndx, nTasks,
							reps, 0);

    // Determine order for printing lines of output
    unsigned ord[npat];
    orderpat(npat, ord, pat);

    printf("#%14s %10s\n", "SitePat", "Prob");
    char        buff[100];
    for(j = 0; j < npat; ++j) {
        char        buff2[100];
        snprintf(buff2, sizeof(buff2), "%s",
                 patLbl(sizeof(buff), buff, pat[ord[j]],
                        GPTree_getLblNdxPtr(gptree)));
        printf("%15s %10.7lf\n", buff2, prob[ord[j]]);
    }

    GPTree_free(gptree);

    return 0;
}
