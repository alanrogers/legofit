/**
 * @file gptree.c
 * @author Alan R. Rogers
 * @brief Methods for simulating gene genealogies within a given tree
 * of populations, and allowing populations to mix and also to split.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "gene.h"
#include "gptree.h"
#include "lblndx.h"
#include "parse.h"
#include "param.h"
#include "parstore.h"
#include "ptrptrmap.h"
#include "sampndx.h"
#include <errno.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#include <unistd.h>

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
    PopNode    *rootPop;        // root of population tree
    Gene       *rootGene;       // root of gene tree
    Bounds      bnd;            // legal range of twoN parameters and time pars
    ParStore   *parstore;       // All parameters
    LblNdx      lblndx;         // Index of sample labels
    SampNdx     sndx;           // Index of sample pointers into PopNode objs
};

int         GPTree_equals(const GPTree *lhs, const GPTree *rhs);

/// Initialize vector x. If ndx==0, simply copy the parameter vector
/// from the GPTree object. Otherwise, randomize the GPTree first.
/// This ensures that differential evolution starts with a set of
/// points, one of which is the same as the values in the input
/// file. This allows you to improve on existing estimates without
/// starting from scratch each time.
void GPTree_initStateVec(void *gpt, int ndx, int n, double x[n],
                         gsl_rng *rng){
    if(ndx == 0)
        GPTree_getParams(gpt, n, x);
    else {
        GPTree *g2 = GPTree_dup(gpt);
        GPTree_randomize(g2, rng);
        GPTree_getParams(g2, n, x);
        GPTree_free(g2);
    }
}

/// Print a description of parameters.
void GPTree_printParStore(void * vself, FILE * fp) {
    GPTree *self = vself;
    ParStore_print(self->parstore, fp);
}

/// Print a description of free parameters.
void GPTree_printParStoreFree(void * vself, FILE * fp) {
    GPTree *self = vself;
    ParStore_printFree(self->parstore, fp);
    ParStore_printConstrained(self->parstore, fp);
}

/// Return pointer to name of i'th free parameter
const char *GPTree_getNameFree(void * self, int i) {
    return ParStore_getNameFree(((GPTree *) self)->parstore, i);
}

/// Randomly perturb all free parameters while maintaining inequality
/// constraints.
void GPTree_randomize(void * vself, gsl_rng * rng) {
    GPTree *self = vself;
    for(int i=0; i < ParStore_nPar(self->parstore); ++i) {
        Param *par = ParStore_getParam(self->parstore, i);
        if(!Param_isFree(par))
            continue;
        double orig = Param_getValue(par);
        double trial = Param_getTrialValue(par, rng);
        if(Param_setValue(par, trial)) {
            fprintf(stderr,"%s:%d: illegal trial value: %lg\n",
                    __FILE__,__LINE__,trial);
            exit(EXIT_FAILURE);
        }
        PopNode_update(self->rootPop, self->parstore);

        // Bisect to satisfy inequality constraints.
        while( !GPTree_feasible(self, 0) ) {
            trial = orig + 0.5*(trial - orig);
            if(Param_setValue(par, trial)) {
                fprintf(stderr,"%s:%d: illegal trial value: %lg\n",
                        __FILE__,__LINE__,trial);
                exit(EXIT_FAILURE);
            }
            PopNode_update(self->rootPop, self->parstore);
        }
    }
}

/// Set free parameters from an array.
/// @param[in] n number of parameters in array, which should equal the
/// number of free parameters in the GPTree.
/// @param[in] x array of parameter values.
/// @return 0 on success, 1 if values of x violate boundary constraints.
int GPTree_setParams(void * vself, int n, double x[n]) {
    GPTree *self = vself;
    assert(n == ParStore_nFree(self->parstore));
    int status = ParStore_setFreeParams(self->parstore, n, x);
    if(status)
        return status;
    PopNode_update(self->rootPop, self->parstore);
    return 0;
}

/// Copy free parameters from GPTree into an array
/// @param[out] n number of parameters in array, which should equal the
/// number of free parameters in the GPTree.
/// @param[out] x array into which parameters will be copied
void GPTree_getParams(void * self, int n, double x[n]) {
    ParStore_getFreeParams(((GPTree*)self)->parstore, n, x);
}

/// Return number of free parameters
int GPTree_nFree(const void * self) {
    return ParStore_nFree(((const GPTree *) self)->parstore);
}

/// Use coalescent simulation to estimate the mean branch length
/// associated with each site pattern.
/// @param self GPTree object
/// @param[out] branchtab BranchTab object, which will tabulate branch
/// lengths from this simulations.
/// @param[inout] rng GSL random number generator
/// @param[in] nreps number of replicate gene trees to simulate
/// @param[in] doSing if doSing is non-zero, singleton site patterns
/// will be tabulated.
void GPTree_brlen(void * vself, BranchTab * branchtab, gsl_rng * rng,
                  unsigned long nreps, int doSing) {
    GPTree *self = vself;
    unsigned long rep;
    for(rep = 0; rep < nreps; ++rep) {
        PopNode_clear(self->rootPop);   // remove old samples

        // Put samples into the gene tree. This allocates memory for
        // each Gene in the sample and puts pointers to them into the
        // PopNodes that are controlled by the SampNdx. The Gene
        // objects aren't owned by SampNdx or PopNode. They will
        // eventually be freed by a call to Gene_free, which
        // recursively frees the root and all descendants.
        for(unsigned i=0; i < SampNdx_size(&self->sndx); ++i) {
            PopNode *node = (PopNode *) SampNdx_get(&self->sndx, i);
            PopNode_newSample(node, i);
        }

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

    BranchTab_divideBy(branchtab, (double) nreps);
}

/// GPTree constructor
void *GPTree_new(const char *fname, Bounds bnd) {
    GPTree     *self = malloc(sizeof(GPTree));
    CHECKMEM(self);

    memset(self, 0, sizeof(GPTree));
    self->bnd = bnd;
    self->parstore = NULL;
    LblNdx_init(&self->lblndx);
    SampNdx_init(&self->sndx);
    FILE       *fp = efopen(fname, "r");
    self->nseg = countSegments(fp);
    rewind(fp);

    PtrPair pp = mktree(fp, &self->sndx, &self->lblndx, &self->bnd);

    self->rootPop = pp.a;
    self->parstore = pp.b;

    fclose(fp);
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
void GPTree_free(void * vself) {
    GPTree *self = vself;
    Gene_free(self->rootGene);
    self->rootGene = NULL;
    PopNode_clear(self->rootPop);
    PopNode_free(self->rootPop);
    ParStore_free(self->parstore);
    free(self);
}

/// Duplicate a GPTree object
void *GPTree_dup(const void * vold) {
    const GPTree *old = vold;
    assert(old);
    if(!GPTree_feasible(old, 1)) {
        pthread_mutex_lock(&outputLock);
        fflush(stdout);
        dostacktrace(__FILE__, __LINE__, stderr);
        fprintf(stderr, "%s:%s:%d: old tree isn't feasible\n", __FILE__,
                __func__,__LINE__);
        ParStore_print(old->parstore, stderr);
        pthread_mutex_unlock(&outputLock);
        exit(EXIT_FAILURE);
    }
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

    GPTree *new = memdup(old, sizeof(GPTree));
    new->parstore = ParStore_dup(old->parstore);
    CHECKMEM(new->parstore);

    assert(old->nseg == new->nseg);

    new->sndx = old->sndx;

    // A hashmap that maps old PopNode pointers to new ones. It is
    // populated by PopNode_dup and used by SampNdx_remapPtrs.
    PtrPtrMap *ppm = PtrPtrMap_new(2 * old->nseg);
    CHECKMEM(ppm);

    new->rootPop = PopNode_dup(old->rootPop, ppm);
    CHECKMEM(new->rootPop);

    SampNdx_remapPtrs(&new->sndx, ppm);

    PtrPtrMap_free(ppm);
    
    assert(PopNode_isClear(new->rootPop));

    GPTree_sanityCheck(new, __FILE__, __LINE__);
    assert(GPTree_equals(old, new));
    if(!GPTree_feasible(new, 1)) {
        pthread_mutex_lock(&outputLock);
        fflush(stdout);
        dostacktrace(__FILE__, __LINE__, stderr);
        fprintf(stderr, "%s:%s:%d: new tree isn't feasible\n", __FILE__,
                __func__,__LINE__);
        ParStore_print(new->parstore, stderr);
        pthread_mutex_unlock(&outputLock);
        exit(EXIT_FAILURE);
    }
    return (void *) new;
}

void GPTree_sanityCheck(void * vself, const char *file, int line) {
#ifndef NDEBUG
    GPTree *self = vself;
    REQUIRE(self->nseg > 0, file, line);
    Bounds_sanityCheck(&self->bnd, file, line);
    ParStore_sanityCheck(self->parstore, file, line);
    LblNdx_sanityCheck(&self->lblndx, file, line);
#endif
}

/// Return 1 if two GPTree objects are equal, 0 if they differ.  Abort
/// with an error if the GPTree pointers are different but one or more
/// of the internal pointers is the same.  Does not access rootPop or
/// rootGene. Used for testing GPTree_dup.
int GPTree_equals(const GPTree * lhs, const GPTree * rhs) {
    if(lhs == rhs)
        return 0;
#ifndef NDEBUG    
    if(lhs->rootPop == rhs->rootPop) {
        fprintf(stderr,
                "%s:%s:%d: two GPTree objects share a rootPop pointer\n",
                __FILE__, __func__, __LINE__);
        exit(EXIT_FAILURE);
    }
    if(lhs->parstore == rhs->parstore) {
        fprintf(stderr,
                "%s:%s:%d: two GPTree objects share a ParStore pointer\n",
                __FILE__, __func__, __LINE__);
        exit(EXIT_FAILURE);
    }
#endif    
    if(!PopNode_equals(lhs->rootPop, rhs->rootPop))
        return 0;
    if(!Bounds_equals(&lhs->bnd, &rhs->bnd))
        return 0;
    if(!ParStore_equals(lhs->parstore, rhs->parstore)) {
        return 0;
    }
    if(!LblNdx_equals(&lhs->lblndx, &rhs->lblndx))
        return 0;
    if(!SampNdx_equals(&lhs->sndx, &rhs->sndx))
        return 0;
    return 1;
}

/// Get the LblNdx object from a GPTree
LblNdx GPTree_getLblNdx(void * vself) {
    return ((GPTree *) vself)->lblndx;
}

/// Are parameters within the feasible region?
int GPTree_feasible(const void * vself, int verbose) {
    const GPTree *self = vself;
    return PopNode_feasible(self->rootPop, self->bnd, verbose);
}

unsigned GPTree_nSamples(void *vself) {
    GPTree *self = vself;

    return SampNdx_size(&self->sndx);
}

#ifdef TEST

#  include <string.h>
#  include <assert.h>
#  include <time.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

#  include "network.h"
#  include <assert.h>
#  include <unistd.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

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
            fprintf(stderr, "usage: xgptree [-v]\n");
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

    Network_init(SIM);

    Bounds      bnd = {
        .lo_twoN = 0.0,
        .hi_twoN = 1e7,
        .lo_t = 0.0,
        .hi_t = INFINITY
    };
    GPTree     *g = GPTree_new(fname, bnd);
    GPTree     *g2 = GPTree_dup(g);
    assert(GPTree_equals(g, g2));

    GPTree_randomize(g2, rng);
    assert( !GPTree_equals(g, g2) );
    gsl_rng_free(rng);
    rng = NULL;

    if(verbose) {
        fprintf(stderr,"Before randomization:\n");
        GPTree_printParStore(g, stderr);
        fprintf(stderr,"After randomization:\n");
        GPTree_printParStore(g2, stderr);
    }

    const LblNdx lblndx = GPTree_getLblNdx(g);
    if(verbose)
        LblNdx_print(&lblndx, stdout);

    GPTree_free(g);
    GPTree_free(g2);

    unlink(fname);
    unitTstResult("GPTree", "OK");
    return 0;
}
#endif
