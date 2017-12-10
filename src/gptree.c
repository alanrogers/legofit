/**
 * @file gptree.c
 * @author Alan R. Rogers
 * @brief Methods for simulating gene genealogies within a given tree
 * of populations, and allowing populations to mix and also to split.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "gptree.h"
#include "gene.h"
#include "lblndx.h"
#include "parse.h"
#include "parstore.h"
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include <pthread.h>
extern pthread_mutex_t outputLock;

/// GPTree stands for Gene-Population tree. It represents a network
/// of populations, which can split to form daughter populations or
/// exchange genes at various points in time. The population tree
/// is static--it does not evolve, but is specified by parsing an
/// external file. Within the population tree, there is a gene
/// tree. The gene tree is dynamic: it can be generated repeatedly via
/// coalescent simulation.
struct GPTree {
    int         nseg;           // number of segments in population tree.
    PopNode    *pnv;            // array of nseg PopNode objects
    PopNode    *rootPop;        // root of population tree
    Gene       *rootGene;       // root of gene tree
    Bounds      bnd;            // legal range of twoN parameters and time parameters
    ParStore   *parstore;       // Fixed and free parameters
    LblNdx      lblndx;         // Index of sample labels
    SampNdx     sndx;           // Index of sample pointers into PopNode objects.
};

/// Print a description of parameters.
void GPTree_printParStore(GPTree * self, FILE * fp) {
    if(ParStore_constrain(self->parstore))
        fprintf(fp,"Warning: free parameters violate constraints\n");
    ParStore_print(self->parstore, fp);
}

/// Print a description of free parameters.
void GPTree_printParStoreFree(GPTree * self, FILE * fp) {
    ParStore_printFree(self->parstore, fp);
    ParStore_printConstrained(self->parstore, fp);
}

/// Randomly perturb all free parameters in the population tree while
/// maintaining inequality constraints.
void GPTree_randomize(GPTree * self, gsl_rng * rng) {
    PopNode_randomize(self->rootPop, self->bnd, self->parstore, rng);
}

/// Set free parameters from an array.
/// @param[in] n number of parameters in array, which should equal the
/// number of free parameters in the GPTree.
/// @param[in] x array of parameter values.
int GPTree_setParams(GPTree * self, int n, double x[n]) {
    assert(n == ParStore_nFree(self->parstore));
    ParStore_setFreeParams(self->parstore, n, x);
    return ParStore_constrain(self->parstore);
}

/// Copy free parameters from GPTree into an array
/// @param[out] n number of parameters in array, which should equal the
/// number of free parameters in the GPTree.
/// @param[out] x array into which parameters will be copied
void GPTree_getParams(GPTree * self, int n, double x[n]) {
    assert(n == ParStore_nFree(self->parstore));
    ParStore_getFreeParams(self->parstore, n, x);
}

/// Return number of free parameters
int GPTree_nFree(const GPTree * self) {
    return ParStore_nFree(self->parstore);
}

/// Build a gene tree by coalescent simulation and then tabulate
/// branch lengths associated with each site pattern.
/// @param self GPTree object
/// @param[out] branchtab BranchTab object, which will tabulate branch
/// lengths from this (and other) simulations.
/// @param[inout] rng GSL random number generator
/// @param[in] nreps number of replicate gene trees to simulate
/// @param[in] doSing if doSing is non-zero, singleton site patterns
/// will be tabulated.
void GPTree_simulate(GPTree * self, BranchTab * branchtab, gsl_rng * rng,
                     unsigned long nreps, int doSing) {
    unsigned long rep;
    if(ParStore_constrain(self->parstore)) {
        fprintf(stderr,"%s:%d: free parameters violate constraints\n",
                __FILE__,__LINE__);
    }
    for(rep = 0; rep < nreps; ++rep) {
        PopNode_clear(self->rootPop);   // remove old samples
        SampNdx_populateTree(&(self->sndx));    // add new samples
        PopNode_gaussian(self->rootPop, self->bnd, self->parstore, rng);

        // coalescent simulation generates gene genealogy within
        // population tree.
        self->rootGene = PopNode_coalesce(self->rootPop, rng);
        assert(self->rootGene);

        // Traverse gene tree, accumulating branch lengths in bins
        // that correspond to site patterns.
        Gene_tabulate(self->rootGene, branchtab, doSing);

        // Free gene genealogy but not population tree.
        Gene_free(self->rootGene);
        self->rootGene = NULL;
    }
}

/// GPTree constructor
GPTree     *GPTree_new(const char *fname, Bounds bnd) {
    GPTree     *self = malloc(sizeof(GPTree));
    CHECKMEM(self);

    memset(self, 0, sizeof(GPTree));
    self->bnd = bnd;
    self->parstore = ParStore_new();
    LblNdx_init(&self->lblndx);
    SampNdx_init(&self->sndx);
    FILE       *fp = efopen(fname, "r");
    self->nseg = countSegments(fp);
    rewind(fp);
    self->pnv = malloc(self->nseg * sizeof(self->pnv[0]));
    CHECKMEM(self->pnv);

    NodeStore  *ns = NodeStore_new(self->nseg, self->pnv);
    CHECKMEM(ns);

    self->rootPop = mktree(fp, &self->sndx, &self->lblndx, self->parstore,
                           &self->bnd, ns);

    fclose(fp);
    NodeStore_free(ns);
    GPTree_sanityCheck(self, __FILE__, __LINE__);
    if(!GPTree_feasible(self, 1)) {
        fprintf(stderr,
                "%s:%s:%d: file \"%s\" describes an infeasible tree.\n",
                __FILE__, __func__, __LINE__, fname);
        GPTree_printParStore(self, stderr);
        exit(EXIT_FAILURE);
    }
    return self;
}

/// GPTree destructor.
void GPTree_free(GPTree * self) {
    Gene_free(self->rootGene);
    self->rootGene = NULL;
    PopNode_clear(self->rootPop);
    self->rootPop = NULL;
    free(self->pnv);
    ParStore_free(self->parstore);
    free(self);
}

/// Duplicate a GPTree object
GPTree     *GPTree_dup(const GPTree * old) {
    assert(old);
    if(ParStore_constrain(old->parstore)) {
        fprintf(stderr,"%s:%d: free parameters violate constraints\n",
                __FILE__,__LINE__);
    }
    assert(GPTree_feasible(old, 0));
    if(old->rootGene != NULL) {
        fprintf(stderr, "%s:%s:%d: old->rootGene must be NULL on entry\n",
                __FILE__, __func__, __LINE__);
        exit(EXIT_FAILURE);
    }
    if(!PopNode_isClear(old->rootPop)) {
        fprintf(stderr, "%s:%s:%d: clear GPTree of samples before call"
                " to GPTree_dup\n", __FILE__, __func__, __LINE__);
        exit(EXIT_FAILURE);
    }

    GPTree     *new = memdup(old, sizeof(GPTree));
    new->parstore = ParStore_dup(old->parstore);
    new->pnv = memdup(old->pnv, old->nseg * sizeof(PopNode));

    assert(old->nseg == new->nseg);
    CHECKMEM(new->parstore);
    CHECKMEM(new->pnv);

    new->sndx = old->sndx;

    /*
     * Adjust the pointers so they refer to the memory allocated in
     * "new" rather than that in "old".  dpar is the absolute offset
     * between new and old for parameter pointers.  dpop is the
     * analogous offset for PopNode pointers. spar and spop are the
     * signs of these offsets. Everything has to be cast to size_t,
     * because we are not doing pointer arithmetic in the usual
     * sense. Ordinarily, ptr+3 means ptr + 3*sizeof(*ptr). We want it
     * to mean ptr+3*sizeof(char).
     */
    int         i, spar, spop;
    size_t      dpar, dpop;

    // Calculate offsets and signs
    if(new->parstore > old->parstore) {
        dpar = ((size_t) new->parstore) - ((size_t) old->parstore);
        spar = 1;
    } else {
        dpar = ((size_t) old->parstore) - ((size_t) new->parstore);
        spar = -1;
    }
    if(new->pnv > old->pnv) {
        dpop = ((size_t) new->pnv) - ((size_t) old->pnv);
        spop = 1;
    } else {
        dpop = ((size_t) old->pnv) - ((size_t) new->pnv);
        spop = -1;
    }

    SHIFT_PTR(new->rootPop, dpop, spop);
    for(i = 0; i < old->nseg; ++i) {
        PopNode_shiftParamPtrs(&new->pnv[i], dpar, spar);
        PopNode_shiftPopNodePtrs(&new->pnv[i], dpop, spop);
    }
    SampNdx_shiftPtrs(&new->sndx, dpop, spop);
    assert(SampNdx_ptrsLegal(&new->sndx, new->pnv, new->pnv + new->nseg));

    GPTree_sanityCheck(new, __FILE__, __LINE__);
    assert(GPTree_equals(old, new));
    if(!GPTree_feasible(new, 0)) {
        pthread_mutex_lock(&outputLock);
        fflush(stdout);
        dostacktrace(__FILE__, __LINE__, stderr);
        fprintf(stderr, "%s:%d: new tree isn't feasible\n", __FILE__,
                __LINE__);
        pthread_mutex_unlock(&outputLock);
        exit(EXIT_FAILURE);
    }
    if(ParStore_constrain(new->parstore)) {
        fprintf(stderr,"%s:%d: free parameters violate constraints\n",
                __FILE__,__LINE__);
    }
    assert(GPTree_feasible(new, 1));
    return new;
}

void GPTree_sanityCheck(GPTree * self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self->nseg > 0, file, line);
    REQUIRE(self->pnv != NULL, file, line);
    REQUIRE(self->rootPop >= self->pnv, file, line);
    REQUIRE(self->rootPop < self->pnv + self->nseg, file, line);
    Bounds_sanityCheck(&self->bnd, file, line);
    ParStore_sanityCheck(self->parstore, file, line);
    LblNdx_sanityCheck(&self->lblndx, file, line);
#endif
}

/// Return 1 if two GPTree objects are equal, 0 if they differ.  Abort
/// with an error if the GPTree pointers are different but one or more
/// of the internal pointers is the same.  Does not access rootPop or
/// rootGene.
int GPTree_equals(const GPTree * lhs, const GPTree * rhs) {
    if(lhs == rhs)
        return 0;
    if(lhs->pnv == rhs->pnv)
        eprintf("%s:%s:%d: two GPTree objects share a PopNode pointer\n",
                __FILE__, __func__, __LINE__);
    if(lhs->parstore == rhs->parstore)
        eprintf("%s:%s:%d: two GPTree objects share a ParStore pointer\n",
                __FILE__, __func__, __LINE__);
    if(!Bounds_equals(&lhs->bnd, &rhs->bnd))
        return 0;
    if(!ParStore_equals(lhs->parstore, rhs->parstore)) {
        fprintf(stderr,"%s:%d: !ParStore_equals\n",__FILE__,__LINE__);
        return 0;
    }
    if(!LblNdx_equals(&lhs->lblndx, &rhs->lblndx))
        return 0;
    if(!SampNdx_equals(&lhs->sndx, &rhs->sndx))
        return 0;
    return 1;
}

/// Get the LblNdx object from a GPTree
LblNdx GPTree_getLblNdx(GPTree * self) {
    return self->lblndx;
}

/// Return pointer to array of lower bounds of free parameters
double     *GPTree_loBounds(GPTree * self) {
    return ParStore_loBounds(self->parstore);
}

/// Return pointer to array of upper bounds of free parameters
double     *GPTree_upBounds(GPTree * self) {
    return ParStore_upBounds(self->parstore);
}

/// Return number of samples.
unsigned GPTree_nsamples(GPTree * self) {
    return SampNdx_size(&self->sndx);
}

/// Are parameters within the feasible region?
int GPTree_feasible(const GPTree * self, int verbose) {
    if(ParStore_constrain(self->parstore))
        return 0;
    return PopNode_feasible(self->rootPop, self->bnd, verbose);
}

#ifdef TEST

#  include <string.h>
#  include <assert.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

#  include <assert.h>
#  include <unistd.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

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
    "time free   Tc=1\n"
    "time free   Tab=3\n"
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
            fprintf(stderr, "usage: xgptree [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    const char *fname = "mktree-tmp.lgo";
    FILE       *fp = fopen(fname, "w");
    fputs(tstInput, fp);
    fclose(fp);

    Bounds      bnd = {
        .lo_twoN = 0.0,
        .hi_twoN = 1e7,
        .lo_t = 0.0,
        .hi_t = HUGE_VAL
    };
    GPTree     *g = GPTree_new(fname, bnd);
    GPTree     *g2 = GPTree_dup(g);
    assert(GPTree_equals(g, g2));

    const LblNdx lblndx = GPTree_getLblNdx(g);
    if(verbose)
        LblNdx_print(&lblndx, stdout);

    GPTree_free(g);
    GPTree_free(g2);

    unlink(fname);
    unitTstResult("GPTree", "untested");
    return 0;
}
#endif
