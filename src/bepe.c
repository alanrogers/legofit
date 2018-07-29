/**
@file bepe.c
@page bepe
@author Daniel R. Tabin and Alan R. Rogers
@brief Bootstrap estimate of prediction error.

# `bepe`: calculate the bootstrap estimate of prediction error

This program, like it's sibling @ref clic "clic", provides a tool for
selecting among models that differ in complexity. It implements the
"bootstrap estimate of predictive error", which is described in
section 17.6 of *An introduction to the bootstrap*, (Efron and
Tibshirani, 1993).  This method provides a solution to the problem of
@ref modsel "overfitting".

Bepe does this by using bootstrap replicates as a proxy for samples
from the underlying (and usually unknown) statistical
distribution. Its value estimates the mean squared error between
predicted site pattern frequencies and those of unobserved samples
from the same statistical distribution. The best model is the one for
which bepe reports the smallest value.

    usage: bepe <realdat> <bdat1> <bdat2> ... -L <real.legofit> 
     <b1.legofit> <b2.legofit> ...

     where realdat is the real data, each "bdat" file is the data
     for one bootstrap replicate, and each "b#.legofit" file is the
     legofit output from the corresponding bootstrap replicate
     Must include realdat file and at least 2 bootstrap replicates.

    Options:
       -h or --help
       print this message

In typical usage, one would type something like

    bepe realdat.txt boot*.txt -L real.legofit boot*.legofit

This usage assumes that your computer's shell or command interpreter sorts
the files globbed by `boot*.txt` and `boot*.legofit` in a consistent order,
so that the i'th .legofit file is the output produced from the i'th
bootstrap data file.

@copyright Copyright (c) 2018, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "strdblqueue.h"
#include "misc.h"
#include <strings.h>

// Function prototypes
void usage(void);
static void checkConsistency(const char *fname1, const char *fname2,
                             StrDblQueue *q1, StrDblQueue *q2);

//vars
const char *usageMsg =
    "usage: bepe <realdat> <bdat1> <bdat2> ... -L <real.legofit>\n"
    " <b1.legofit> <b2.legofit> ...\n" "\n"
    " where realdat is the real data, each \"bdat\" file is the data\n"
    " for one bootstrap replicate, and each \"b#.legofit\" file is the\n"
    " legofit output from the corresponding bootstrap replicate\n"
    " Must include realdat file and at least 2 bootstrap replicates.\n" "\n"
    "Options:\n" "   -h or --help\n" "   print this message\n";
void usage(void) {
    fputs(usageMsg, stderr);
    exit(EXIT_FAILURE);
} 

static void checkConsistency(const char *fname1, const char *fname2,
                             StrDblQueue *q1, StrDblQueue *q2) {
    if(StrDblQueue_compare(q1, q2)) {
        fprintf(stderr, "%s:%d: inconsistent site patterns in"
                " files %s and %s\n", __FILE__, __LINE__,
                fname1, fname2);
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char **argv) {

    // Command line arguments specify file names
    int i, j;
    int nfiles=0, nLegoFiles=0;
    int gotDashL = 0;

    for(i = 1; i < argc; i++) {
        if(argv[i][0] == '-') {
            if(strcmp(argv[i], "-L") == 0) {
                gotDashL = 1;
                continue;
            }else
                usage();
        }
        if(gotDashL)
            ++nLegoFiles;
        else
            ++nfiles;
    }
    if(nfiles != nLegoFiles) {
        fprintf(stderr, "%s:%d\n"
                " Inconsistent number of files!"
                " %d data files and %d legofit files\n",
                __FILE__, __LINE__, nfiles, nLegoFiles);
        usage();
    }

    if(nfiles < 3)
        usage();

    const char *datafname[nfiles], *legofname[nfiles];
    gotDashL=0;
    j=0;
    for(i = 1; i < argc; i++) {
        if(argv[i][0] == '-') {
            if(strcmp(argv[i], "-L") == 0) {
                gotDashL = 1;
                assert(j==nfiles);
                j = 0;
                continue;
            }else
                usage();
        }
        if(gotDashL)
            legofname[j++] = argv[i];
        else
            datafname[j++] = argv[i];
    }
    assert(j==nfiles);

    // Read data and bootstrap files into an arrays of queues
    StrDblQueue *data_queue[nfiles];
    StrDblQueue *lego_queue[nfiles];
    for(i = 0; i < nfiles; ++i) {
        lego_queue[i] = StrDblQueue_parseSitPat(legofname[i]);
        data_queue[i] = StrDblQueue_parseSitPat(datafname[i]);
        StrDblQueue_normalize(lego_queue[i]);
        StrDblQueue_normalize(data_queue[i]);
        if(i==0) {
            checkConsistency(datafname[0], legofname[0],
                             data_queue[0], lego_queue[0]);
            continue;
        }
        checkConsistency(legofname[0], legofname[i],
                         lego_queue[0], lego_queue[i]);
        checkConsistency(datafname[0], legofname[i],
                         data_queue[0], data_queue[i]);
    }

    // i is the file against which all others are being compared.
    // j indexes the other files.
    for(i=0; i<nfiles; ++i) {
        double msd = 0.0, bias=0.0, bepe;
        for(j=0; j < nfiles; ++j) {
            if(j==i)
                continue;
            msd += StrDblQueue_msd(data_queue[i], lego_queue[j]);
            bias += StrDblQueue_msd(data_queue[j], lego_queue[j]);
        }
        bepe = (msd+bias)/(nfiles-1);
        printf("%lg \t#BEPE based on %s\n", bepe, mybasename(datafname[i]));
    }

    // Echo input files
    printf("# %2s %15s %15s\n", "i", "datfile[i]", "legofile[i]");
    for(i = 0; i < nfiles; ++i) {
        printf("# %2d %15s %15s\n", i,
                mybasename(datafname[i]),
                mybasename(legofname[i]));
    }

    // Echo command line
    printf("# cmd:");
    for(i = 0; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');
    fflush(stdout);

    return 0;
}
