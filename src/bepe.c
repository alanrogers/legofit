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

Bepe does this by using bootstrap or simulation replicates as a proxy
for samples from the underlying (and usually unknown) statistical
distribution. Its value estimates the mean squared error between
predicted site pattern frequencies and those of unobserved samples
from the same statistical distribution. The best model is the one for
which bepe reports the smallest value.

    usage: bepe <realdat> <bdat1> <bdat2> ... -L <real.legofit>
     <b1.legofit> <b2.legofit> ...

     where realdat is the real data, each "bdat" file is the data
     for one bootstrap or simulation replicate, and each "b#.legofit"
     file is the legofit output from the corresponding bootstrap
     replicate Must include realdat file and at least 2 bootstrap
     replicates.

    Options:
       -h or --help   : print this message

In typical usage, one would type something like

    bepe realdat.txt boot*.txt -L real.legofit boot*.legofit

This usage assumes that your computer's shell or command interpreter sorts
the files globbed by `boot*.txt` and `boot*.legofit` in a consistent order,
so that the i'th .legofit file is the output produced from the i'th
bootstrap data file.

The lines of data in `bepe`'s output are arranged according to the
order in which data files are listed on the command line. It is
important that this order be consistent, if multiple `.bepe` files are
to be compared using @ref booma "booma". (Otherwise, `booma` will
abort with an error.) Consistency problems can arise because of
differences in locale settings on different machines. The example
command in the preceding paragraph may generate file names in
different orders on different machines. To ensure a consistent order
on Unix-like operating systems (linux and osx), set the LC_ALL
environment parameter first. Using the bash shell:

    export LC_ALL=C
    bepe realdat.txt boot*.txt -L real.legofit boot*.legofit

This ensures that file names will be listed in the same order,
regardless of the locale setting of the local machine.

The first few lines of `bepe`'s output look like this:

    ################################################
    # bepe: bootstrap estimate of predictive error #
    #                 version 1.80                 #
    ################################################

    # Program was compiled: Feb  8 2019 15:26:39
    # Program was run: Fri Feb  8 15:31:50 2019

    #          bepe        DataFile     LegofitFile
    1.908861378e-07        sim0.opf    a2-0.legofit
    2.495970126e-07        sim1.opf    a2-1.legofit
    2.627267214e-07       sim10.opf   a2-10.legofit
     1.71297571e-07       sim11.opf   a2-11.legofit
    1.917256551e-07       sim12.opf   a2-12.legofit

Many of the lines begin with the comment character (#) so that they
can be removed easily in downstream processing. Uncommented lines
refer to particular data sets, whose names are given in columns
2-3. The order of these lines reflects the order in which files are
listed on the command line. If you are analyzing real data plus
bootstrap replicates, it is useful to list the real data file first on
the command line, so that its bebe value will appear first in the
output. Column 1 (bepe) gives the mean squared difference between the
site pattern frequencies in this data set, and the predictions
obtained from all the other data sets. Bepe also includes a bias
correction, which is described by Efron and Tibshirani. The predicted
values are taken from the various `.legofit` files.

@copyright Copyright (c) 2018, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "strdblqueue.h"
#include "misc.h"
#include <strings.h>
#include <time.h>

// Function prototypes
void usage(void);
static void checkConsistency(const char *fname1, const char *fname2,
                             StrDblQueue *q1, StrDblQueue *q2);

//vars
const char *usageMsg =
    "usage: bepe [options] <realdat> <bdat1> <bdat2> ... -L <real.legofit>\n"
    " <b1.legofit> <b2.legofit> ...\n" "\n"
    " where realdat is the real data, each \"bdat\" file is the data\n"
    " for one bootstrap or simulation replicate, and each \"b#.legofit\"\n"
    " file is the legofit output from the corresponding bootstrap replicate\n"
    " Must include realdat file and at least 2 bootstrap replicates.\n" "\n"
    "Options:\n"
    "   -h or --help   : print this message\n";

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

    time_t currtime = time(NULL);

    hdr("bepe: bootstrap estimate of predictive error");
#if defined(__DATE__) && defined(__TIME__)
    printf("# Program was compiled: %s %s\n", __DATE__, __TIME__);
#endif
    printf("# Program was run: %s", ctime(&currtime));

    int i, j;
    int nfiles=0, nLegoFiles=0;
    int gotDashL = 0;

    for(i = 1; i < argc; i++) {
        if(argv[i][0] == '-') {
            if(strcmp(argv[i], "-L") == 0) {
                gotDashL = 1;
                continue;
            }else {
                fprintf(stderr,"unknown flag argument: %s\n", argv[i]);
                usage();
            }
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

    if(nfiles < 3) {
        fprintf(stderr,"nfiles=%d; need at least 3\n", nfiles);
        usage();
    }

    const char *datafname[nfiles], *legofname[nfiles];
    gotDashL=0;
    j=0;
    for(i = 1; i < argc; i++) {
        if(argv[i][0] == '-') {
            if(strcmp(argv[i], "-L") == 0) {
                gotDashL = 1;
                assert(j==nfiles);
                j = 0;
            }
            continue;
        }
        if(gotDashL)
            legofname[j++] = argv[i];
        else
            datafname[j++] = argv[i];
    }
    assert(j==nfiles);

    putchar('\n');

    // Read data and legofit files into an arrays of queues
    StrDblQueue *data_queue[nfiles], *lego_queue[nfiles];
    for(i = 0; i < nfiles; ++i) {
        lego_queue[i] = StrDblQueue_parseSitePat(legofname[i]);
        data_queue[i] = StrDblQueue_parseSitePat(datafname[i]);
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
    printf("#%14s %15s %15s\n", "bepe", "DataFile", "LegofitFile");
    for(i=0; i<nfiles; ++i) {
        double msd = 0.0, bias=0.0, bepe;
        for(j=0; j < nfiles; ++j) {
            if(j==i)
                continue;
            msd += StrDblQueue_msd(data_queue[i], lego_queue[j]);
            bias += StrDblQueue_msd(data_queue[j], lego_queue[j]);
        }
        bepe = (msd+bias)/(nfiles-1);
        printf("%15.10lg %15s %15s\n", bepe,
               mybasename(datafname[i]),
               mybasename(legofname[i]));
    }
    putchar('\n');

    return 0;
}
