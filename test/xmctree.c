/**
 * @file xmctree.c
 * @author Alan R. Rogers
 * @brief Test mctree.c.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "lblndx.h"
#include "mctree.h"
#include "misc.h"
#include "parstore.h"
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <unistd.h>

#ifdef NDEBUG
#  error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include "network.h"
#include <assert.h>
#include <unistd.h>

#ifdef NDEBUG
#  error "Unit tests must be compiled without -DNDEBUG flag"
#endif

//      a-------|
//              |ab--|
//      b--|b2--|    |
//         |         |abc--
//         |c--------|
//
//  t = 0  1    3    5.5     inf
const char *tstInput =
    " # this is a comment\n"
    "time fixed  T0=0\n"
    "time free   x = 2\n"
    "time free   Tc=1\n"
    "time constrained Tab=x - Tc\n"
    "time free   Tabc=5.5\n"
    "twoN free   twoNa=100\n"
    "twoN fixed  twoNb=123\n"
    "twoN free   twoNc=213.4\n"
    "twoN fixed  twoNb2=32.1\n"
    "twoN free   twoNab=222\n"
    "twoN fixed  twoNabc=1.2e2\n"
    "mixFrac free Mc=0.02\n"
    "segment a   t=T0     twoN=twoNa    samples=1\n"
    "segment b   t=T0     twoN=twoNb    samples=1\n"
    "segment c   t=Tc     twoN=twoNc    samples=1\n"
    "segment b2  t=Tc     twoN=twoNb2\n"
    "segment ab  t=Tab    twoN=twoNab\n"
    "segment abc t=Tabc   twoN=twoNabc\n"
    "mix    b  from b2 + Mc * c\n"
    "derive a  from ab\n"
    "derive b2 from ab\n" "derive ab from abc\n" "derive c  from abc\n";
int main(int argc, char **argv) {

    int         verbose = 0;

    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xmctree [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    long unsigned seed = time(NULL);
    seed ^= getpid();
    gsl_rng_set(rng, seed);

    const char *fname = "mktree-tmp.lgo";
    FILE       *fp = fopen(fname, "w");
    fputs(tstInput, fp);
    fclose(fp);

    Network_init(DETERMINISTIC);

    Bounds      bnd = {
        .lo_twoN = 0.0,
        .hi_twoN = 1e7,
        .lo_t = 0.0,
        .hi_t = INFINITY
    };
    MCTree     *g = MCTree_new(fname, bnd);

    MCTree     *g2 = MCTree_dup(g);
    assert(MCTree_equals(g, g2));

    MCTree_randomize(g2, rng);
    assert( !MCTree_equals(g, g2) );
    gsl_rng_free(rng);
    rng = NULL;

    if(verbose) {
        fprintf(stderr,"Before randomization:\n");
        MCTree_printParStore(g, stderr);
        fprintf(stderr,"After randomization:\n");
        MCTree_printParStore(g2, stderr);
    }

    const LblNdx lblndx = MCTree_getLblNdx(g);
    if(verbose)
        LblNdx_print(&lblndx, stdout);

    MCTree_free(g);
    MCTree_free(g2);

    unlink(fname);
    unitTstResult("MCTree", "OK");
    return 0;
}
