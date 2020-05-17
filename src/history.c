/**
 * @file history.c
 * @author Alan R. Rogers
 * @brief Generic interface for GPtree and MCTree.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "history.h"
#include "gptree.h"
#include "mctree.h"
#include <gsl/gsl_rng.h>

struct History {
    int method;  // either MATCOAL or SIM
    GPTree *gptree;
    MCTree *mctree;
};

void History_randomize(History * self, gsl_rng * rng);

/// Initialize vector x.
void History_initStateVec(History *self, int ndx, int n, double x[n],
                         gsl_rng *rng){
    switch(self->method) {
    case SIM:
        GPTree_initStateVec(self->gptree, ndx, n, x, rng);
        break;
    case MATCOAL:
        MCTree_initStateVec(self->mctree, ndx, n, x, rng);
        break;
    default:
        fprintf(stderr,"%s:%d: bad method\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
}

/// Print a description of parameters.
void History_printParStore(History * self, FILE * fp) {
    switch(self->method) {
    case SIM:
        GPTree_printParStore(self->gptree, fp);
        break;
    case MATCOAL:
        MCTree_printParStore(self->mctree, fp);
        break;
    default:
        fprintf(stderr,"%s:%d: bad method\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
}

/// Print a description of free parameters.
void History_printParStoreFree(History * self, FILE * fp) {
    switch(self->method) {
    case SIM:
        GPTree_printParStoreFree(self->gptree, fp);
        break;
    case MATCOAL:
        MCTree_printParStoreFree(self->mctree, fp);
        break;
    default:
        fprintf(stderr,"%s:%d: bad method\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
}

/// Return pointer to name of i'th free parameter
const char *History_getNameFree(History * self, int i) {
    switch(self->method) {
    case SIM:
        return GPTree_getNameFree(self->gptree, i);
    case MATCOAL:
        return MCTree_getNameFree(self->mctree, i);
    default:
        fprintf(stderr,"%s:%d: bad method\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
}

/// Randomly perturb all free parameters while maintaining inequality
/// constraints.
void History_randomize(History * self, gsl_rng * rng) {
    switch(self->method) {
    case SIM:
        GPTree_randomize(self->gptree, rng);
        break;
    case MATCOAL:
        MCTree_randomize(self->mctree, rng);
        break;
    default:
        fprintf(stderr,"%s:%d: bad method\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
}


/// Set free parameters from an array.
/// @param[in] n number of parameters in array, which should equal the
/// number of free parameters in the GPTree.
/// @param[in] x array of parameter values.
/// @return 0 on success, 1 if values of x violate boundary constraints.
int History_setParams(History * self, int n, double x[n]) {
    switch(self->method) {
    case SIM:
        return GPTree_setParams(self->gptree, n, x);
    case MATCOAL:
        return MCTree_setParams(self->mctree, n, x);
    default:
        fprintf(stderr,"%s:%d: bad method\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
}

/// Copy free parameters from GPTree into an array
/// @param[out] n number of parameters in array, which should equal the
/// number of free parameters in the GPTree.
/// @param[out] x array into which parameters will be copied
void History_getParams(History * self, int n, double x[n]) {
    switch(self->method) {
    case SIM:
        GPTree_getParams(self->gptree, n, x);
        break;
    case MATCOAL:
        MCTree_getParams(self->mctree, n, x);
        break;
    default:
        fprintf(stderr,"%s:%d: bad method\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
}

/// Return number of free parameters
int History_nFree(const History * self) {
    switch(self->method) {
    case SIM:
        return GPTree_nFree(self->gptree);
    case MATCOAL:
        return MCTree_nFree(self->mctree);
    default:
        fprintf(stderr,"%s:%d: bad method\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
}

/// Estimate the probability of each site pattern.
/// @param self GPTree object
/// @param[out] branchtab BranchTab object, which will tabulate branch
/// lengths from this simulations.
/// @param[inout] rng GSL random number generator
/// @param[in] nreps number of replicate gene trees to simulate
/// @param[in] doSing if doSing is non-zero, singleton site patterns
/// will be tabulated.
void History_patprob(History * self, BranchTab * branchtab, gsl_rng * rng,
                     unsigned long nreps, int doSing) {
    switch(self->method) {
    case SIM:
        GPTree_patprob(self->gptree, branchtab, rng, nreps, doSing);
        break;
    case MATCOAL:
        MCTree_patprob(self->mctree, branchtab, rng, nreps, doSing);
        break;
    default:
        fprintf(stderr,"%s:%d: bad method\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
}

/// constructor
History *History_new(const char *fname, Bounds bnd, int method) {
    History *self = malloc(sizeof(History));
    CHECKMEM(self);

    self->method = method;
    switch(self->method) {
    case SIM:
        self->gptree = GPTree_new(fname, bnd);
        self->mctree = NULL;
        break;
    case MATCOAL:
        self->gptree = NULL;
        self->mctree = MCTree_new(fname, bnd);
        break;
    default:
        fprintf(stderr,"%s:%d: bad method\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    return self;
}

/// destructor.
void History_free(History * self) {
    if(self->gptree)
        GPTree_free(self->gptree);
    if(self->mctree)
        MCTree_free(self->mctree);
    free(self);
}

/// Duplicate
History *History_dup(const History * old) {
    if(old == NULL)
        return NULL;
    History *new = malloc(sizeof(History));
    CHECKMEM(new);

    new->method = old->method;
    switch(new->method) {
    case SIM:
        new->gptree = GPTree_dup(old->gptree);
        new->mctree = NULL;
        break;
    case MATCOAL:
        new->gptree = NULL;
        new->mctree = MCTree_dup(old->mctree);
        break;
    default:
        fprintf(stderr,"%s:%d: bad method\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    return new;
}

void History_sanityCheck(History * self, const char *file, int line) {
#ifndef NDEBUG
    switch(self->method) {
    case SIM:
        GPTree_sanityCheck(self->gptree, file, line);
        REQUIRE(self->mctree == NULL, file, line);
        break;
    case MATCOAL:
        REQUIRE(self->gptree == NULL, file, line);
        MCTree_sanityCheck(self->mctree, file, line);
        break;
    default:
        fprintf(stderr,"%s:%d: bad method\n",file,line);
        exit(EXIT_FAILURE);
    }
#endif
}

/// Get the LblNdx object
LblNdx History_getLblNdx(History * self) {
    switch(self->method) {
    case SIM:
        return GPTree_getLblNdx(self->gptree);
    case MATCOAL:
        return MCTree_getLblNdx(self->mctree);
    default:
        fprintf(stderr,"%s:%d: bad method\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
}

/// Are parameters within the feasible region?
int History_feasible(const History * self, int verbose) {
    switch(self->method) {
    case SIM:
        return GPTree_feasible(self->gptree, verbose);
    case MATCOAL:
        return MCTree_feasible(self->mctree, verbose);
    default:
        fprintf(stderr,"%s:%d: bad method\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
}

#ifdef TEST

#  include <string.h>
#  include <assert.h>
#  include <time.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

#  include <assert.h>
#  include <unistd.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

int History_equals(History *a, History *b);
int History_equals(History *a, History *b) {
    switch(self->method) {
    case SIM:
        return GPTree_equals(a->gptree, b->gptree);
    case MATCOAL:
        return MCTree_equals(a->mctree, b->mctree)
    default:
        fprintf(stderr,"%s:%d: bad method\n",__FILE__,__LINE__);
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
            fprintf(stderr, "usage: xhistory [-v]\n");
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
    int method = SIM;
    History    *g = History_new(fname, bnd, method);
    History    *g2 = History_dup(g);
    assert(History_equals(g, g2));

    History_randomize(g2, rng);
    assert( !History_equals(g, g2) );
    gsl_rng_free(rng);
    rng = NULL;

    if(verbose) {
        fprintf(stderr,"Before randomization:\n");
        History_printParStore(g, stderr);
        fprintf(stderr,"After randomization:\n");
        History_printParStore(g2, stderr);
    }

    const LblNdx lblndx = History_getLblNdx(g);
    if(verbose)
        LblNdx_print(&lblndx, stdout);

    History_free(g);
    History_free(g2);

    unlink(fname);
    unitTstResult("History", "OK");
    return 0;
}
#endif
