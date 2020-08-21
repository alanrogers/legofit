/**
 * @file network.c
 * @author Alan R. Rogers
 * @brief Models a network of populations
 *
 * Legofit supports two models of the network of populations. One of
 * these (PopNode and GPTree) calculate probabilities by coalescent
 * simulation. The other (Segment and MCTree) does the same
 * calculation via a deterministic algorithm. This file implements a
 * generic interface to these two models.
 * 
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "typedefs.h"
#include "gptree.h"
#include "popnode.h"
#include "segment.h"
#include "mctree.h"
#include <gsl/gsl_rng.h>

// External variables defined here
#define EXTERN
#include "network.h"
#undef EXTERN

struct Network {
    enum NetworkType method;  // either DETERMINISTIC or STOCHASTIC
    GPTree *gptree;
    MCTree *mctree;
};

/// Initialize model of network
void Network_init(enum NetworkType type) {
    switch(type) {
    case STOCHASTIC:
        Network_dup = GPTree_dup;
        Network_feasible = GPTree_feasible;
        Network_free = GPTree_free;
        Network_getLblNdx = GPTree_getLblNdx;
        Network_getNameFree = GPTree_getNameFree;
        Network_getParams = GPTree_getParams;
        Network_new = GPTree_new;
        Network_nFree = GPTree_nFree;
        Network_patprob = GPTree_patprob;
        Network_printParStore = GPTree_printParStore;
        Network_randomize = GPTree_randomize;
        Network_sanityCheck = GPTree_sanityCheck;
        Network_setParams = GPTree_setParams;
        Network_initStateVec = GPTree_initStateVec;
        Node_new = PopNode_new;
        Node_addChild = PopNode_addChild;
        Node_mix = PopNode_mix;
        Node_root = PopNode_root;
        Node_print = PopNode_print;
        break;
    case DETERMINISTIC:
        Network_dup = MCTree_dup;
        Network_feasible = MCTree_feasible;
        Network_free = MCTree_free;
        Network_getLblNdx = MCTree_getLblNdx;
        Network_getNameFree = MCTree_getNameFree;
        Network_getParams = MCTree_getParams;
        Network_new = MCTree_new;
        Network_nFree = MCTree_nFree;
        Network_patprob = MCTree_patprob;
        Network_printParStore = MCTree_printParStore;
        Network_randomize = MCTree_randomize;
        Network_sanityCheck = MCTree_sanityCheck;
        Network_setParams = MCTree_setParams;
        Network_initStateVec = MCTree_initStateVec;
        Node_new = Segment_new;
        Node_addChild = Segment_addChild;
        Node_mix = Segment_mix;
        Node_root = Segment_root;
        Node_print = Segment_print;
        break;
    default:
        fprintf(stderr,"%s:%d: unknown type (%d)\n",
                __FILE__,__LINE__, type);
        exit(EXIT_FAILURE);
    }
}

#ifdef TEST

#  include <string.h>
#  include <assert.h>
#  include <time.h>
#  include <unistd.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

int Network_equals(Network *a, Network *b);

int Network_equals(Network *a, Network *b) {
    switch(self->method) {
    case STOCHASTIC:
        return GPTree_equals(a->gptree, b->gptree);
    case DETERMINISTIC:
        return MCTree_equals(a->mctree, b->mctree)
    default:
        fprintf(stderr,"%s:%d: bad method (%d)\n",
                __FILE__,__LINE__, self->method);
        exit(EXIT_FAILURE);
    }
}


//      a-------|
//              |ab--|
//      b--|bb--|    |
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
    "twoN fixed  twoNbb=32.1\n"
    "twoN free   twoNab=222\n"
    "twoN fixed  twoNabc=1.2e2\n"
    "mixFrac free Mc=0.02\n"
    "segment a   t=T0     twoN=twoNa    samples=1\n"
    "segment b   t=T0     twoN=twoNb    samples=1\n"
    "segment c   t=Tc     twoN=twoNc    samples=1\n"
    "segment bb  t=Tc     twoN=twoNbb\n"
    "segment ab  t=Tab    twoN=twoNab\n"
    "segment abc t=Tabc   twoN=twoNabc\n"
    "mix    b  from bb + Mc * c\n"
    "derive a  from ab\n"
    "derive bb from ab\n" "derive ab from abc\n" "derive c  from abc\n";
int main(int argc, char **argv) {

    int         verbose = 0;

    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xnetwork [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, (unsigned long) time(NULL));

    const char *fname = "mktree-tmp.lgo";
    FILE       *fp = fopen(fname, "w");
    fputs(tstInput, fp);
    fclose(fp);

    Bounds      bnd = {
        .lo_twoN = 0.0,
        .hi_twoN = 1e7,
        .lo_t = 0.0,
        .hi_t = INFINITY
    };
    int method = STOCHASTIC;
    Network    *g = Network_new(fname, bnd, method);
    Network    *g2 = Network_dup(g);
    assert(Network_equals(g, g2));

    Network_randomize(g2, rng);
    assert( !Network_equals(g, g2) );
    gsl_rng_free(rng);
    rng = NULL;

    if(verbose) {
        fprintf(stderr,"Before randomization:\n");
        Network_printParStore(g, stderr);
        fprintf(stderr,"After randomization:\n");
        Network_printParStore(g2, stderr);
    }

    const LblNdx lblndx = Network_getLblNdx(g);
    if(verbose)
        LblNdx_print(&lblndx, stdout);

    Network_free(g);
    Network_free(g2);

    unlink(fname);
    unitTstResult("Network", "OK");
    return 0;
}
#endif
