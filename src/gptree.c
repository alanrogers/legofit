/**
 * @file gptree.c
 * @brief Methods for simulating gene genealogies within a given tree
 * of populations, and allowing populations to mix and also to split.
 */

#include "gptree.h"
#include "gene.h"
#include "lblndx.h"
#include "parse.h"
#include "parstore.h"
#include <string.h>

struct GPTree {
    int nseg; // number of segments in population tree.
    PopNode *pnv; // array of length nseg
    PopNode *rootPop;
    Gene *rootGene;
    int nNodes;
    Bounds bnd;
    ParStore *parstore;
    LblNdx lblndx;
    SampNdx sndx;
};

GPTree *GPTree_new(const char *fname, Bounds bnd) {
    GPTree *self = malloc(sizeof(GPTree));
    CHECKMEM(self);

    memset(self, 0, sizeof(GPTree));
    self->bnd = bnd;
    self->parstore = ParStore_new();
    LblNdx_init(&self->lblndx);
    SampNdx_init(&self->sndx);
    FILE *fp = fopen(fname, "r");
    if(fp == NULL)
        eprintf("%s:%s:%d: can't open file \"%s\".\n",
                __FILE__, __func__, __LINE__, fname);
    self->nseg = countSegments(fp);
    rewind(fp);
    self->pnv = malloc(self->nseg * sizeof(self->pnv[0]));
    CHECKMEM(self->pnv);

    NodeStore *ns = NodeStore_new(self->nseg, self->pnv);
    CHECKMEM(ns);

    self->rootPop = mktree(fp, &self->sndx, &self->lblndx, self->parstore,
                           &self->bnd, ns);

    fclose(fp);
    NodeStore_free(ns);
    return self;
}

/// Destructor.
void GPTree_free(GPTree *self) {
    Gene_free(self->rootGene);
    PopNode_clear(self->rootPop);
    free(self->pnv);
    ParStore_free(self->parstore);
    free(self);
}

/// Duplicate a GPTree object
GPTree *GPTree_dup(GPTree *old) {
    Gene_free(old->rootGene);
    old->rootGene = NULL;
    PopNode_clear(old->rootPop);

    GPTree *new   = memdup(old, sizeof(GPTree));
    new->parstore = ParStore_dup(old->parstore);
    new->pnv      = memdup(old->pnv, old->nseg * sizeof(PopNode));
    assert(old->nseg == new->nseg);

    /*
     * Adjust the pointers within each PopNode object so they refer to
     * the memory allocated in "new" rather than that in "old".
     * dpar is the offset between new and old for parameter pointers.
     * dpop is the analogous offset for PopNode pointers. Everything
     * has to be cast to size_t, because we are not doing pointer
     * arithmetic in the usual sense. Ordinarily, ptr+3 means
     * ptr + 3*sizeof(*ptr). We want it to mean ptr+3*sizeof(char).
     */
    int i;
    size_t dpar = ((size_t) new->parstore) - ((size_t) old->parstore);
    size_t dpop = ((size_t) new->pnv) - ((size_t) old->pnv);
    for(i=0; i < old->nseg; ++i) {
        PopNode_shiftParamPtrs(&new->pnv[i], dpar);
        PopNode_shiftPopNodePtrs(&new->pnv[i], dpop);
    }

    return new;
}

#ifdef TEST

#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

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
    "twoN free   2Na=100\n"
    "twoN fixed  2Nb=123\n"
    "twoN free   2Nc=213.4\n"
    "twoN fixed  2Nbb=32.1\n"
    "twoN free   2Nab=222\n"
    "twoN fixed  2Nabc=1.2e2\n"
    "mixFrac free Mc=0.02\n"
    "segment a   t=T0     twoN=2Na    samples=1\n"
    "segment b   t=T0     twoN=2Nb    samples=1\n"
    "segment c   t=Tc     twoN=2Nc    samples=1\n"
    "segment bb  t=Tc     twoN=2Nbb\n"
    "segment ab  t=Tab    twoN=2Nab\n"
    "segment abc t=Tabc   twoN=2Nabc\n"
    "mix    b  from bb + Mc * c\n"
    "derive a  from ab\n"
    "derive bb from ab\n"
    "derive ab from abc\n"
    "derive c  from abc\n";
int main(int argc, char **argv) {

    int verbose=0;

    if(argc > 1) {
        if(argc!=2 || 0!=strcmp(argv[1], "-v")) {
            fprintf(stderr,"usage: xgptree [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    const char *tstFname = "mktree-tmp.lgo";
    FILE       *fp = fopen(tstFname, "w");
    fputs(tstInput, fp);
    fclose(fp);
    fp = fopen(tstFname, "r");
    if(fp == NULL) {
        fprintf(stderr, "%s:%d: Can't open file \"%s\"\n",
                __FILE__, __LINE__, tstFname);
        exit(1);
    }

    fclose(fp);
    unitTstResult("GPTree", "untested");
    return 0;
}
#endif
