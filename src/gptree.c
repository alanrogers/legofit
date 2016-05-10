/**
 * @file gptree.c
 * @brief Methods for simulating gene genealogies within a given tree
 * of populations, and allowing populations to mix and also to split.
 */

#include "gptree.h"
#include "gene.h"
#include "misc.h"
#include "binary.h"
#include "branchtab.h"
#include "parstore.h"
#include "parse.h"
#include "lblndx.h"
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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

struct NodeStore {
    int nused, len;
    PopNode *v; // not locally owned
};


/// Set everything to zero.
void SampNdx_init(SampNdx * self) {
    memset(self, 0, sizeof(*self));
}

/// Add samples for a single population. Should be called once for
/// each sampled population.
void SampNdx_addSamples(SampNdx * self, unsigned nsamples,
						PopNode * pnode) {
    unsigned    i;
    if(self->n + nsamples >= MAXSAMP)
        eprintf("%s:%s:%d: too many samples\n", __FILE__, __func__, __LINE__);
    for(i = 0; i < nsamples; ++i) {
        self->node[self->n] = pnode;
        self->n += 1;
    }
}

/// Put samples into the gene tree. Should be done at the start of
/// each simulation.
void SampNdx_populateTree(SampNdx * self) {
    unsigned    i;
    for(i = 0; i < self->n; ++i)
        PopNode_newGene(self->node[i], i);
}

unsigned SampNdx_size(SampNdx * self) {
    return self->n;
}

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

// Add increment INC to pointer PTR. Units are sizeof(char)
// rather than the size of the object to which PTR refers.
#define INCR_PTR(PTR,INC) do{                         \
        (PTR) = ((size_t) (PTR)) + ((size_t) (INC));  \
    }while(0);

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
    int i, j;
    size_t dpar = ((size_t) new->parstore) - ((size_t) old->parstore);
    size_t dpop = ((size_t) new->pnv) - ((size_t) old->pnv);
    for(i=0; i < old->nseg; ++i) {
        INCR_PTR(new->pnv[i].twoN,  dpar);
        INCR_PTR(new->pnv[i].start, dpar);
        INCR_PTR(new->pnv[i].end,   dpar);
        for(j = 0; j < new->nparents; ++j)
            INCR_PTR(new->pnv[i].parent[j], dpop);
        for(j = 0; j < new->nchildren; ++j)
            INCR_PTR(new->pnv[i].child[j],  dpop);
    }

    return new;
}

NodeStore *NodeStore_new(int len, PopNode *v) {
    NodeStore *self = malloc(sizeof(NodeStore));
    CHECKMEM(self);

    self->nused = 0;
    self->len = len;
    self->v = v;
    return self;
}

void NodeStore_free(NodeStore *self) {
    // Does not free self->v
    free(self);
}

PopNode *NodeStore_alloc(NodeStore *self) {
    if(self->nused >= self->len)
        eprintf("%s:%s:%d: Ran out of PopNode objects.\n",
                __FILE__, __func__, __LINE__);
    return &self->v[self->nused++];
}

#ifdef TEST

#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {

    int verbose=0;

    if(argc > 1) {
        if(argc!=2 || 0!=strcmp(argv[1], "-v")) {
            fprintf(stderr,"usage: xgptree [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    tipId_t id1 = 1;
    tipId_t id2 = 2;

    int nseg = 10;
    PopNode v[nseg];
    NodeStore *ns = NodeStore_new(nseg, v);
    CHECKMEM(ns);

    PopNode *node = NodeStore_alloc(ns);
    assert(node == v);
    node = NodeStore_alloc(ns);
    assert(node == &v[1]);
    node = NodeStore_alloc(ns);
    assert(node == &v[2]);

    NodeStore_free(ns);
	unitTstResult("NodeStore", "OK");
    
    ns = NodeStore_new(nseg, v);
    CHECKMEM(ns);

    double twoN0 = 1.0, start0= 0.0;
    PopNode *p0 = PopNode_new(&twoN0, &start0, ns);
    assert(p0->twoN == &twoN0);
    assert(p0->start == &start0);
    assert(p0->end == NULL);
    assert(p0->mix == NULL);
    assert(p0->nsamples == 0);
    assert(p0->nchildren == 0);
    assert(p0->child[0] == NULL);
    assert(p0->child[1] == NULL);
    assert(p0->parent[0] == NULL);
    assert(p0->parent[1] == NULL);

	if(verbose) 
		PopNode_printShallow(p0, stdout);

    double twoN1 = 100.0, start1= 123.0;
    PopNode *p1 = PopNode_new(&twoN1, &start1, ns);
    assert(p1->twoN == &twoN1);
    assert(p1->start == &start1);
    assert(p1->end == NULL);
    assert(p1->mix == NULL);
    assert(p1->nsamples == 0);
    assert(p1->nchildren == 0);
    assert(p1->child[0] == NULL);
    assert(p1->child[1] == NULL);
    assert(p1->parent[0] == NULL);
    assert(p1->parent[1] == NULL);

    Gene *g1 = Gene_new(id1);
    Gene *g2 = Gene_new(id2);
    PopNode_addSample(p1, g1);
    PopNode_addSample(p1, g2);
    assert(p1->nsamples == 2);

    PopNode_addChild(p1, p0);
    assert(p1->nchildren == 1);
    assert(p0->nparents == 1);
    assert(p1->child[0] == p0);
    assert(p0->parent[0] == p1);

	if(verbose) 
		PopNode_printShallow(p1, stdout);

    unitTstResult("PopNode", "untested");

    SampNdx     sndx = {.n = 3 };
    assert(sndx.n == 3);
    assert(SampNdx_size(&sndx) == 3);

    SampNdx_init(&sndx);
    assert(SampNdx_size(&sndx) == 0);

	double twoN = 100.0;
	double start = 20.0;
    PopNode    *pnode = PopNode_new(&twoN, &start, ns);
    SampNdx_addSamples(&sndx, 1, pnode);
    SampNdx_addSamples(&sndx, 2, pnode);

    assert(3 == SampNdx_size(&sndx));
    SampNdx_populateTree(&sndx);
    assert(3 == PopNode_nsamples(pnode));
    NodeStore_free(ns);

	unitTstResult("SampNdx", "OK");

    return 0;
}
#endif
