/**
 * @file gptree.c
 * @author Alan R. Rogers
 * @brief Methods for simulating gene genealogies within a given tree
 * of populations, and allowing populations to mix and also to split.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#define UNUSED(x) (void)(x)

#include "branchtab.h"
#include "lblndx.h"
#include "matcoal.h"
#include "mctree.h"
#include "param.h"
#include "parse.h"
#include "parstore.h"
#include "ptrptrmap.h"
#include "sampndx.h"
#include "segment.h"
#include <string.h>
#include <pthread.h>

extern pthread_mutex_t outputLock;

/// MCTree stands for Gene-Population tree. It represents a network
/// of populations, which can split to form daughter populations or
/// exchange genes at various points in time. The population tree
/// is static--it does not evolve, but is specified by parsing an
/// external file. Within the population tree, there is a gene
/// tree. The gene tree is dynamic: it can be generated repeatedly via
/// coalescent simulation.
struct MCTree {
    int         nseg;           // number of segments in population tree.
    Segment    *rootPop;        // root of population tree
    Bounds      bnd;            // legal range of twoN parameters and time pars
    ParStore   *parstore;       // All parameters
    LblNdx      lblndx;         // Index of sample labels
    SampNdx     sndx;           // Index of sample pointers into Segment objs
};

/// MCTree constructor
void *MCTree_new(const char *fname, Bounds bnd) {
    MCTree     *self = malloc(sizeof(MCTree));
    CHECKMEM(self);

    memset(self, 0, sizeof(MCTree));
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

    // Initialize MatCoal.
    unsigned nsamples = SampNdx_size(&self->sndx);
    int mcsamples = MatCoal_nSamples();
    if(mcsamples == 0) {
        MatCoal_initExterns(nsamples);
    }else if(mcsamples < nsamples) {
        fprintf(stderr,"%s:%d: MatCoal has been initialized with %d samples,"
                " but file %s specifies a larger number, %u.\n",
                __FILE__,__LINE__, mcsamples, fname, nsamples);
        exit(EXIT_FAILURE);
    }

    fclose(fp);
    MCTree_sanityCheck(self, __FILE__, __LINE__);
    if(!MCTree_feasible(self, 1)) {
        fprintf(stderr,
                "%s:%s:%d: file \"%s\" describes an infeasible tree.\n",
                __FILE__, __func__, __LINE__, fname);
        MCTree_printParStore(self, stderr);
        exit(EXIT_FAILURE);
    }
    return self;
}

/// MCTree destructor.
void MCTree_free(void * vself) {
    MCTree *self = vself;
    Segment_free(self->rootPop);
    ParStore_free(self->parstore);
    free(self);
}

void MCTree_sanityCheck(void * vself, const char *file, int line) {
#ifndef NDEBUG
    MCTree *self = vself;
    REQUIRE(self->nseg > 0, file, line);
    Bounds_sanityCheck(&self->bnd, file, line);
    ParStore_sanityCheck(self->parstore, file, line);
    LblNdx_sanityCheck(&self->lblndx, file, line);
    Segment_sanityCheck(self->rootPop, file, line);
#endif
}

/// Are parameters within the feasible region?
int MCTree_feasible(const void * vself, int verbose) {
    const MCTree *self = vself;
    return Segment_feasible(self->rootPop, self->bnd, verbose);
}

/// Print a description of parameters.
void MCTree_printParStore(void * vself, FILE * fp) {
    MCTree *self = vself;
    ParStore_print(self->parstore, fp);
}

/// Get the LblNdx object from a MCTree
LblNdx MCTree_getLblNdx(void * vself) {
    return ((MCTree *) vself)->lblndx;
}

/// Duplicate a MCTree object
void *MCTree_dup(const void * vold) {
    const MCTree *old = vold;
    assert(old);
    if(!MCTree_feasible(old, 1)) {
        pthread_mutex_lock(&outputLock);
        fflush(stdout);
        dostacktrace(__FILE__, __LINE__, stderr);
        fprintf(stderr, "%s:%s:%d: old tree isn't feasible\n", __FILE__,
                __func__,__LINE__);
        ParStore_print(old->parstore, stderr);
        pthread_mutex_unlock(&outputLock);
        exit(EXIT_FAILURE);
    }
    if(!Segment_isClear(old->rootPop)) {
        fprintf(stderr, "%s:%s:%d: clear MCTree of samples before call"
                " to MCTree_dup\n", __FILE__, __func__, __LINE__);
        exit(EXIT_FAILURE);
    }

    MCTree *new = memdup(old, sizeof(MCTree));
    new->parstore = ParStore_dup(old->parstore);
    CHECKMEM(new->parstore);

    assert(old->nseg == new->nseg);

    new->sndx = old->sndx;

    // A hashmap that maps old Segment pointers to new ones. It is
    // populated by Segment_dup and used by SampNdx_remapPtrs.
    PtrPtrMap *ppm = PtrPtrMap_new();
    CHECKMEM(ppm);

    new->rootPop = Segment_dup(old->rootPop, ppm);
    CHECKMEM(new->rootPop);

    SampNdx_remapPtrs(&new->sndx, ppm);

    PtrPtrMap_free(ppm);
    
    assert(Segment_isClear(new->rootPop));

    MCTree_sanityCheck(new, __FILE__, __LINE__);
    assert(MCTree_equals(old, new));
    if(!MCTree_feasible(new, 1)) {
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

/// Return 1 if two MCTree objects are equal, 0 if they differ.  Abort
/// with an error if the MCTree pointers are different but one or more
/// of the internal pointers is the same.  Does not access rootPop or
/// rootGene. Used for testing MCTree_dup.
int MCTree_equals(const MCTree * lhs, const MCTree * rhs) {
    if(lhs == rhs)
        return 0;
#ifndef NDEBUG    
    if(lhs->rootPop == rhs->rootPop) {
        fprintf(stderr,
                "%s:%s:%d: two MCTree objects share a rootPop pointer\n",
                __FILE__, __func__, __LINE__);
        exit(EXIT_FAILURE);
    }
    if(lhs->parstore == rhs->parstore) {
        fprintf(stderr,
                "%s:%s:%d: two MCTree objects share a ParStore pointer\n",
                __FILE__, __func__, __LINE__);
        exit(EXIT_FAILURE);
    }
#endif    
    if(!Segment_equals(lhs->rootPop, rhs->rootPop))
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

/// Randomly perturb all free parameters while maintaining inequality
/// constraints.
void MCTree_randomize(void * vself, gsl_rng * rng) {
    MCTree *self = vself;
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
        // Bisect to satisfy inequality constraints.
        while( !MCTree_feasible(self, 0) ) {
            trial = orig + 0.5*(trial - orig);
            if(Param_setValue(par, trial)) {
                fprintf(stderr,"%s:%d: illegal trial value: %lg; type=%d\n",
                        __FILE__,__LINE__,trial, par->type);
                exit(EXIT_FAILURE);
            }
        }
    }
}

/// Return pointer to name of i'th free parameter
const char *MCTree_getNameFree(void * self, int i) {
    return ParStore_getNameFree(((MCTree *) self)->parstore, i);
}

/// Copy free parameters from MCTree into an array
/// @param[out] n number of parameters in array, which should equal the
/// number of free parameters in the MCTree.
/// @param[out] x array into which parameters will be copied
void MCTree_getParams(void * self, int n, double x[n]) {
    ParStore_getFreeParams(((MCTree*)self)->parstore, n, x);
}

/// Return number of free parameters
int MCTree_nFree(const void * self) {
    return ParStore_nFree(((const MCTree *) self)->parstore);
}

/// Set free parameters from an array.
/// @param[in] n number of parameters in array, which should equal the
/// number of free parameters in the MCTree.
/// @param[in] x array of parameter values.
/// @return 0 on success, 1 if values of x violate boundary constraints.
int MCTree_setParams(void * vself, int n, double x[n]) {
    MCTree *self = vself;
    assert(n == ParStore_nFree(self->parstore));
    return ParStore_setFreeParams(self->parstore, n, x);
}

/// Initialize vector x. If ndx==0, simply copy the parameter vector
/// from the MCTree object. Otherwise, randomize the MCTree first.
/// This ensures that differential evolution starts with a set of
/// points, one of which is the same as the values in the input
/// file. This allows you to improve on existing estimates without
/// starting from scratch each time.
void MCTree_initStateVec(void *self, int ndx, int n, double x[n],
                         gsl_rng *rng){
    if(ndx == 0)
        MCTree_getParams(self, n, x);
    else {
        MCTree *g2 = MCTree_dup(self);
        MCTree_randomize(g2, rng);
        MCTree_getParams(g2, n, x);
        MCTree_free(g2);
    }
}

/// Calculate the probability of each site pattern.
/// @param self MCTree object
/// @param[out] branchtab BranchTab object, which will tabulate branch
/// lengths from this simulations.
/// @param[inout] rng GSL random number generator (NOT USED)
/// @param[in] nreps number of replicate gene trees to simulate (NOT USED)
/// @param[in] doSing if doSing is non-zero, singleton site patterns
/// will be tabulated.
void MCTree_patprob(void * vself, BranchTab * branchtab, gsl_rng * rng,
                    unsigned long nreps, int doSing) {
    UNUSED(rng);
    UNUSED(nreps);
    MCTree *self = vself;
    Segment_clear(self->rootPop);   // remove old samples

    // Put samples into the gene tree. This allocates memory for
    // each Gene in the sample and puts pointers to them into the
    // Segments that are controlled by the SampNdx. The Gene
    // objects aren't owned by SampNdx or Segment. They will
    // eventually be freed by a call to Gene_free, which
    // recursively frees the root and all descendants.
    for(unsigned i=0; i < SampNdx_size(&self->sndx); ++i) {
        Segment *node = (Segment *) SampNdx_get(&self->sndx, i);
        Segment_newSample(node, i);
    }

    // Use deterministic algorithm to calculate expected branch lengths
    // for each site pattern.
    int status = Segment_coalesce(self->rootPop, doSing, branchtab);
    if(status) {
        fprintf(stderr,"%s:%d: Segment_coalesce returned status %d\n",
                __FILE__, __LINE__, status);
        exit(EXIT_FAILURE);
    }
    //BranchTab_normalize(branchtab);
}

unsigned MCTree_nSamples(void *vself) {
    MCTree *self = vself;

    return SampNdx_size(&self->sndx);
}

