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

const char *tstInput[7] = {
    // 0: good input
    "chr1" TAB "1" TAB "a" TAB "t" TAB "0/0" TAB "0|1" "\n"
    "chr1" TAB "5" TAB "G" TAB "c" TAB "0/1" TAB "0|1" "\n"
    "chr2" TAB "1" TAB "a" TAB "t" TAB "0/0" TAB "0|1" "\n"
    "chr2" TAB "5" TAB "G" TAB "c" TAB "0/1" TAB "0|1" "\n",

    // 1: Missorted chromosomes
    "chr2" TAB "1" TAB "a" TAB "t" TAB "0/0" TAB "0|1" "\n"
    "chr1" TAB "5" TAB "G" TAB "c" TAB "0/1" TAB "0|1" "\n",

    // 2: Duplicate position
    "chr1" TAB "1" TAB "G" TAB "c" TAB "0/1" TAB "0|1" "\n"
    "chr1" TAB "1" TAB "a" TAB "t" TAB "0/0" TAB "0|1" "\n",

    // 3: Missorted positions
    "chr1" TAB "5" TAB "G" TAB "c" TAB "0/1" TAB "0|1" "\n"
    "chr1" TAB "1" TAB "a" TAB "t" TAB "0/0" TAB "0|1" "\n",

    // 4: Missing field
    "chr1" TAB "1" TAB TAB "t" TAB "0/0" TAB "0|1" "\n",

    // 5: Bad genotype
    "chr1" TAB "1" TAB "a" TAB "t" TAB "0/0" TAB "0" "\n",

    // 6: Bad genotype
    "chr1" TAB "1" TAB "a" TAB "t" TAB "0/0" TAB "0/a" "\n",
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

    // input file 0
    ifp = fopen(infname, "w");
    assert(ifp);
    fputs(tstInput[0], ifp);
    fclose(ifp);
    ifp = fopen(infname, "r");
    assert(ifp);
    ofp = fopen(outfname, "w");
    if(ofp == NULL) {
        fprintf(stderr, "%s:%d: can't write file %s\n",
                __FILE__, __LINE__, outfname);
        exit(1);
    }
    status = raf(ifp, ofp, 0);
    assert(status==0);
    fclose(ifp);
    fclose(ofp);
    RAFReader  *rdr = RAFReader_new(outfname);
    assert(rdr);
    assert(0 == RAFReader_next(rdr));
    assert(RAFReader_raf(rdr) == 0.75);

    assert(0 == RAFReader_next(rdr));
    assert(RAFReader_raf(rdr) == 0.5);

    assert(0 == RAFReader_next(rdr));
    assert(RAFReader_raf(rdr) == 0.75);

    assert(0 == RAFReader_next(rdr));
    assert(RAFReader_raf(rdr) == 0.5);

    assert(EOF == RAFReader_next(rdr));
    RAFReader_free(rdr);

    // input file 1
    ifp = fopen(infname, "w");
    assert(ifp);
    fputs(tstInput[1], ifp);
    fclose(ifp);
    ifp = fopen(infname, "r");
    assert(ifp);
    ofp = fopen(outfname, "w");
    if(ofp == NULL) {
        fprintf(stderr, "%s:%d: can't write file %s\n",
                __FILE__, __LINE__, outfname);
        exit(1);
    }
    status = raf(ifp, ofp, 0);
    assert(status==BAD_SORT);
    fclose(ifp);
    fclose(ofp);

    // input file 2
    ifp = fopen(infname, "w");
    assert(ifp);
    fputs(tstInput[2], ifp);
    fclose(ifp);
    ifp = fopen(infname, "r");
    assert(ifp);
    ofp = fopen(outfname, "w");
    if(ofp == NULL) {
        fprintf(stderr, "%s:%d: can't write file %s\n",
                __FILE__, __LINE__, outfname);
        exit(1);
    }
    status = raf(ifp, ofp, 0);
    assert(status==DUPLICATE_NUCPOS);
    fclose(ifp);
    fclose(ofp);

    // input file 3
    ifp = fopen(infname, "w");
    assert(ifp);
    fputs(tstInput[3], ifp);
    fclose(ifp);
    ifp = fopen(infname, "r");
    assert(ifp);
    ofp = fopen(outfname, "w");
    if(ofp == NULL) {
        fprintf(stderr, "%s:%d: can't write file %s\n",
                __FILE__, __LINE__, outfname);
        exit(1);
    }
    status = raf(ifp, ofp, 0);
    assert(status==BAD_SORT);
    fclose(ifp);
    fclose(ofp);

    // input file 4
    ifp = fopen(infname, "w");
    assert(ifp);
    fputs(tstInput[4], ifp);
    fclose(ifp);
    ifp = fopen(infname, "r");
    assert(ifp);
    ofp = fopen(outfname, "w");
    if(ofp == NULL) {
        fprintf(stderr, "%s:%d: can't write file %s\n",
                __FILE__, __LINE__, outfname);
        exit(1);
    }
    status = raf(ifp, ofp, 0);
    if(status!=0) {
        fprintf(stderr,"%s:%d: bad rtn from raf. %d instead of %d\n",
                __FILE__, __LINE__, status, 0);
        unitTstResult("raf", "FAIL");
        exit(1);
    }
    fclose(ifp);
    fclose(ofp);
    
    // input file 5
    ifp = fopen(infname, "w");
    assert(ifp);
    fputs(tstInput[5], ifp);
    fclose(ifp);
    ifp = fopen(infname, "r");
    assert(ifp);
    ofp = fopen(outfname, "w");
    if(ofp == NULL) {
        fprintf(stderr, "%s:%d: can't write file %s\n",
                __FILE__, __LINE__, outfname);
        exit(1);
    }
    status = raf(ifp, ofp, 0);
    if(status!=BAD_GTYPE) {
        fprintf(stderr,"%s:%d: bad rtn from raf. %d instead of %d\n",
                __FILE__, __LINE__, status, BAD_GTYPE);
        unitTstResult("raf", "FAIL");
        exit(1);
    }
    fclose(ifp);
    fclose(ofp);
    
    unitTstResult("raf", "OK");

    return 0;
}
