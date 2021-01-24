/**
   @file xparseopf.c
   @brief Test parseopf.c.

   @copyright Copyright (c) 2020, Alan R. Rogers
   <rogers@anthro.utah.edu>. This file is released under the Internet
   Systems Consortium License, which can be found in file "LICENSE".
*/
#include "parseopf.h"
#include "lblndx.h"
#include "misc.h"
#include "branchtab.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

const char *tstPatProbInput =
    "#SitePat   obs\n"
    "a          1.0\n"
    "b          2.0\n"
    "c          3.0\n"
    "a:b        4.0\n"
    "a:c        5.0\n"
    "b:c        6.0\n"
    "a:b:c      7.0\n";

int main(void) {

    const char *tstPatProbFname = "patprob-tmp.opf";
    FILE *fp = fopen(tstPatProbFname, "w");
    fputs(tstPatProbInput, fp);
    fclose(fp);

    const tipId_t unity = 1u;
    tipId_t a_id = unity;
    tipId_t b_id = unity << 1;
    tipId_t c_id = unity << 2;
    tipId_t ab_id = a_id | b_id;
    tipId_t ac_id = a_id | c_id;
    tipId_t bc_id = b_id | c_id;
    tipId_t abc_id = a_id | b_id | c_id;

    LblNdx lblndx;
    LblNdx_init(&lblndx);
    LblNdx_addSamples(&lblndx, 1, "a");
    LblNdx_addSamples(&lblndx, 1, "b");
    LblNdx_addSamples(&lblndx, 1, "c");

    assert(a_id == LblNdx_getTipId(&lblndx, "a"));
    assert(b_id == LblNdx_getTipId(&lblndx, "b"));
    assert(c_id == LblNdx_getTipId(&lblndx, "c"));
    assert(ab_id == LblNdx_getTipId(&lblndx, "a:b"));
    assert(ac_id == LblNdx_getTipId(&lblndx, "a:c"));
    assert(bc_id == LblNdx_getTipId(&lblndx, "b:c"));
    assert(abc_id == LblNdx_getTipId(&lblndx, "a:b:c"));
    assert(0 == LblNdx_getTipId(&lblndx, "not_there"));

    BranchTab *bt = parseOpf(tstPatProbFname, &lblndx);
    CHECKMEM(bt);

    assert(1.0 == BranchTab_get(bt, a_id));
    assert(2.0 == BranchTab_get(bt, b_id));
    assert(3.0 == BranchTab_get(bt, c_id));
    assert(4.0 == BranchTab_get(bt, ab_id));
    assert(5.0 == BranchTab_get(bt, ac_id));
    assert(6.0 == BranchTab_get(bt, bc_id));
    assert(7.0 == BranchTab_get(bt, abc_id));

    assert(BranchTab_hasSingletons(bt));

    unlink(tstPatProbFname);
    BranchTab_free(bt);

    unitTstResult("parseOpf", "OK");
    return 0;
}    
