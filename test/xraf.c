/**
 * @file xraf.c
 * @author Alan R. Rogers
 * @brief Test rafreader.c.
 * @copyright Copyright (c) 2026, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "raf.h"
#include "rafreader.h"
#include "misc.h"
#include "error.h"
#include <stdio.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#define TAB "\t"

const char *tstInput[6] = {
    // good input
    "chr1" TAB "1" TAB "a" TAB "t" TAB "0/0" TAB "0|1" "\n"
    "chr1" TAB "5" TAB "G" TAB "c" TAB "0/1" TAB "0|1" "\n"
    "chr2" TAB "1" TAB "a" TAB "t" TAB "0/0" TAB "0|1" "\n"
    "chr2" TAB "5" TAB "G" TAB "c" TAB "0/1" TAB "0|1" "\n",

    // Missorted chromosomes
    "chr2" TAB "1" TAB "a" TAB "t" TAB "0/0" TAB "0|1" "\n"
    "chr1" TAB "5" TAB "G" TAB "c" TAB "0/1" TAB "0|1" "\n",

    // Duplicate position
    "chr1" TAB "1" TAB "G" TAB "c" TAB "0/1" TAB "0|1" "\n"
    "chr1" TAB "1" TAB "a" TAB "t" TAB "0/0" TAB "0|1" "\n",

    // Missorted positions
    "chr1" TAB "5" TAB "G" TAB "c" TAB "0/1" TAB "0|1" "\n"
    "chr1" TAB "1" TAB "a" TAB "t" TAB "0/0" TAB "0|1" "\n",

    // Bad genotype
    "chr1" TAB "1" TAB TAB "t" TAB "0/0" TAB "0" "\n",

    // Bad genotype
    "chr1" TAB "1" TAB TAB "t" TAB "0/0" TAB "0/a" "\n",
};

int main(int argc, char **argv) {

    int status;

    if(argc != 1) {
        fprintf(stderr, "usage: xraf [-v]\n");
        exit(1);
    }
    const char *infname = "rafin.tmp";
    const char *outfname = "rafout.tmp";
    FILE *ifp, *ofp;

    // 1st input file
    ifp = fopen(infname, "w");
    assert(ifp);
    fputs(tstInput[0], ifp);
    fclose(ifp);
    ifp = fopen(infname, "r");
    assert(ifp);
    ofp = fopen(outfname, "r");
    assert(ofp);
    status = raf(ifp, ofp);
    assert(status==0);
    fclose(ofp);
    RAFReader  *rdr = RAFReader_new(outfname);
    assert(rdr);
    assert(0 == RAFReader_next(rdr));
    assert(RAFReader_raf(rdr) == 0.25);
        
    assert(0 == RAFReader_next(rdr));
    assert(RAFReader_raf(rdr) == 0.5);
    assert(0 == RAFReader_next(rdr));
    assert(RAFReader_raf(rdr) == 0.25);
    assert(0 == RAFReader_next(rdr));
    assert(RAFReader_raf(rdr) == 0.5);
    assert(EOF == RAFReader_next(rdr));
    RAFReader_free(rdr);

    unitTstResult("raf", "OK");

    return 0;
}
