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

int raf_status(int idata, const char *outfname);

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

    // 4: Missing ref field
    "chr1" TAB "1" TAB TAB "t" TAB "0/0" TAB "0|1" "\n",

    // 5: Bad genotype
    "chr1" TAB "1" TAB "a" TAB "t" TAB "0/0" TAB "0" "\n",

    // 6: Bad genotype
    "chr1" TAB "1" TAB "a" TAB "t" TAB "0/0" TAB "0/a" "\n",

    // 7: Missing genotypes
    "chr1" TAB "1" TAB "a" TAB "t""\n",
};

// Write the dataset indexed by `idata` to disk, process it with
// `raf`, and return the status code.
int raf_status(int idata, const char *outfname) {
    const char *infname = "rafin.tmp";

    FILE *ifp = mustopen(infname, "w", __FILE__, __LINE__);
    fputs(tstInput[idata], ifp);
    fclose(ifp);
    ifp = mustopen(infname, "r", __FILE__, __LINE__);
    FILE *ofp = mustopen(outfname, "w", __FILE__, __LINE__);
    int status = raf(ifp, ofp, 0);
    fclose(ifp);
    fclose(ofp);
    return status;
}

int main(int argc, char **argv) {

    const char *outfname = "rafout.tmp";

    // input file 0
    REQUIRE(0==raf_status(0, outfname), __FILE__, __LINE__);
    RAFReader  *rdr = RAFReader_new(outfname);
    assert(rdr);
    REQUIRE(0 == RAFReader_next(rdr), __FILE__, __LINE__);
    REQUIRE(RAFReader_raf(rdr) == 0.75, __FILE__, __LINE__);

    REQUIRE(0 == RAFReader_next(rdr), __FILE__, __LINE__);
    REQUIRE(RAFReader_raf(rdr) == 0.5, __FILE__, __LINE__);

    REQUIRE(0 == RAFReader_next(rdr), __FILE__, __LINE__);
    REQUIRE(RAFReader_raf(rdr) == 0.75, __FILE__, __LINE__);

    REQUIRE(0 == RAFReader_next(rdr), __FILE__, __LINE__);
    REQUIRE(RAFReader_raf(rdr) == 0.5, __FILE__, __LINE__);

    REQUIRE(EOF == RAFReader_next(rdr), __FILE__, __LINE__);
    RAFReader_free(rdr);

    // input file 1
    REQUIRE(BAD_SORT==raf_status(1, outfname), __FILE__, __LINE__);

    // input file 2
    REQUIRE(DUPLICATE_NUCPOS==raf_status(2, outfname), __FILE__, __LINE__);

    // input file 3
    REQUIRE(BAD_SORT==raf_status(3, outfname), __FILE__, __LINE__);

    // input file 4
    REQUIRE(EMPTY_FIELD==raf_status(4, outfname), __FILE__, __LINE__);
    
    // input file 5
    REQUIRE(BAD_GTYPE==raf_status(5, outfname), __FILE__, __LINE__);
    
    unitTstResult("raf", "OK");

    return 0;
}
