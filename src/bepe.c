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

The first few lines of `bepe`'s output look like this:

    ################################################
    # bepe: bootstrap estimate of predictive error #
    #                 version 1.47                 #
    ################################################
    
    # Program was compiled: Jul 29 2018 08:45:55
    # Program was run: Sun Jul 29 08:53:47 2018
    
    #          bepe DataFile
    3.599051951e-07 mpqva1.txt
    3.937402126e-07 mpqvasing000
     3.65923644e-07 mpqvasing001
    4.205461171e-07 mpqvasing002
    3.841736741e-07 mpqvasing003
    4.103783016e-07 mpqvasing004
    4.357691253e-07 mpqvasing005
    3.996038295e-07 mpqvasing006

Many of the lines begin with the comment character (#) so that they
can be removed easily in downstream processing. Each uncommented line
refers to a particular data set, whose name is given in column
2. Typically, the first file is the real data, and the others are
bootstrap replicates. Column 1 (bepe) gives the mean squared
difference between the site pattern frequencies in this data set, and
the predictions obtained from all the other data sets. (Bepe also
includes a bias correction, which is described by Efron and
Tibshirani.) The predicted values are taken from the various
`.legofit` files.

The next section lists the two sets of input files. Here are a few
illustrative lines:

    #  i         datfile[i]        legofile[i]
    #  0         mpqva1.txt         t3.legofit
    #  1       mpqvasing000    t3boot0.legofit
    #  2       mpqvasing001    t3boot1.legofit
    #  3       mpqvasing002    t3boot2.legofit
    #  4       mpqvasing003    t3boot3.legofit
    #  5       mpqvasing004    t3boot4.legofit

The first column gives the index of the current data set, the second
names the data file containing site pattern frequencie, and the third
names the `legofit` output generated from that data file. The
`legofit` output includes predicted site pattern frequencies, which
are used by `bepe`.

This section of output is provided to help find a particular kind of
mistake. Suppose we had run `bepe` like this:

    bepe mpqva1.txt mpqvasing* -L t3.legofit t3boot*.legofit

This generates erroneous output, which can be detected by examining
the 2nd portion of `bepe`'s output. Here are the first few lines:

    #  i         datfile[i]        legofile[i]
    #  0         mpqva1.txt         t3.legofit
    #  1       mpqvasing000    t3boot0.legofit
    #  2       mpqvasing001   t3boot10.legofit

Note that bootstrap data file `mpqvasing001` is paired with `legofit`
file `t2boot10.legofit`. This is a mistake. This data file should have
been paired with `t2boot1.legofit`. The problem is that the wildcards
characters (\*) expand to file names sorted in lexical order. Because
of the leading zeroes in the data file names, a lexical sort 
puts these names into into numerical order. But there are no leading
zeroes in the `.legofit` filenames, so they sort in a different
order. This generates a mismatch between the data files and the
`.legofit` files.

To avoid such problems, you need to generate both sets of file names
in a consistent order. Here's a `bash` command line that solves the
above problem:

    bepe mpqva1.txt mpqvasing* -L t3.legofit `seq 0 49 | \
    xargs printf "t3boot%d.legofit "`

The last section of `bepe`'s output echoes the command line, which
serves as documentation for this run of `bepe`.

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

    time_t currtime = time(NULL);

    hdr("bepe: bootstrap estimate of predictive error");
#if defined(__DATE__) && defined(__TIME__)
    printf("# Program was compiled: %s %s\n", __DATE__, __TIME__);
#endif
    printf("# Program was run: %s\n", ctime(&currtime));

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
    printf("#%14s %s\n", "bepe", "DataFile");
    for(i=0; i<nfiles; ++i) {
        double msd = 0.0, bias=0.0, bepe;
        for(j=0; j < nfiles; ++j) {
            if(j==i)
                continue;
            msd += StrDblQueue_msd(data_queue[i], lego_queue[j]);
            bias += StrDblQueue_msd(data_queue[j], lego_queue[j]);
        }
        bepe = (msd+bias)/(nfiles-1);
        printf("%15.10lg %s\n", bepe, mybasename(datafname[i]));
    }
    putchar('\n');

    // Echo input files
    printf("# %2s %18s %18s\n", "i", "datfile[i]", "legofile[i]");
    for(i = 0; i < nfiles; ++i) {
        printf("# %2d %18s %18s\n", i,
                mybasename(datafname[i]),
                mybasename(legofname[i]));
    }
    putchar('\n');

    // Echo command line
    printf("# command line:\n#");
    for(i = 0; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');
    fflush(stdout);

    return 0;
}
