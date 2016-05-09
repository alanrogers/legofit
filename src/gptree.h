#ifndef ARR_GPTREE_H
#  define ARR_GPTREE_H

#  include "typedefs.h"
#  include <gsl/gsl_rng.h>

#  define POPNAMESIZE 30
#  define MAXSAMP ((int)(8*sizeof(tipId_t)))

struct SampNdx {
    unsigned    n;              // number of samples
    char        lbl[MAXSAMP][POPNAMESIZE];
    PopNode    *node[MAXSAMP];
};

void        Gene_tabulate(Gene * self, BranchTab * bt);
void        Gene_free(Gene * gene);
PopNode    *PopNode_new(double *twoN, double *start, NodeStore *ns);
void        PopNode_addChild(PopNode * parent, PopNode * child);
void        PopNode_mix(PopNode * child, double *mPtr, PopNode * introgressor,
                        PopNode * native);
void        PopNode_newGene(PopNode * self, unsigned ndx);
void        PopNode_addSample(PopNode * self, Gene * gene);
Gene       *PopNode_coalesce(PopNode * self, gsl_rng * rng);
void        PopNode_free(PopNode * self);
void        PopNode_clear(PopNode * self);
void        PopNode_print(FILE * fp, PopNode * self, int indent);
void        PopNode_printShallow(PopNode * self, FILE * fp);
PopNode    *PopNode_root(PopNode * self);
void        PopNode_sanityFromLeaf(PopNode * self, const char *file,
                                   int line);
int         PopNode_nsamples(PopNode * self);

void        SampNdx_init(SampNdx * self);
void        SampNdx_addSamples(SampNdx * self, unsigned nsamples,
							   PopNode * pnode);
void        SampNdx_populateTree(SampNdx * self);
unsigned    SampNdx_size(SampNdx * self);

GPTree     *GPTree_new(const char *fname, Bounds bnd);
void        GPTree_free(GPTree *self);

NodeStore  *NodeStore_new(int len, PopNode v[len]);
void        NodeStore_free(NodeStore *self);
PopNode    *NodeStore_alloc(NodeStore *self);

#endif
