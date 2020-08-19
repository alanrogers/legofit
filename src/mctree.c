/**
 * @file gptree.c
 * @author Alan R. Rogers
 * @brief Methods for simulating gene genealogies within a given tree
 * of populations, and allowing populations to mix and also to split.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "lblndx.h"
#include "mctree.h"
#include "parse.h"
#include "parstore.h"
#include "sampndx.h"
#include "segment.h"
#include <string.h>

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

