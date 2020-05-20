#ifndef ARR_POPNODE_H
#  define ARR_POPNODE_H

#  include "typedefs.h"
#  include <stdbool.h>
#  include <gsl/gsl_rng.h>

struct PopNode {
    int         nparents, nchildren;
    double      *twoN;           // ptr to current pop size
    double      *start, *end;    // duration of this PopNode
    double      *mix;            // ptr to frac of pop derived from parent[1]
    struct PopNode *parent[2];
    struct PopNode *child[2];

    int         nsamples;
    Gene       *sample[MAXSAMP]; // not locally owned
};

void       *PopNode_new(double *twoN, double *start, NodeStore *ns);
int         PopNode_addChild(void * vparent, void * vchild);
int         PopNode_mix(void * vchild, double *mPtr, void * vintrogressor,
                        void * vnative);
void        PopNode_newGene(PopNode * self, unsigned ndx);
void        PopNode_addSample(PopNode * self, Gene * gene);
Gene       *PopNode_coalesce(PopNode * self, gsl_rng * rng);
void        PopNode_chkDependencies(PopNode * self, ParStore * parstore);
int         PopNode_feasible(const PopNode *self, Bounds bnd, int verbose);
void        PopNode_free(PopNode * self);
void        PopNode_clear(PopNode * self);
int         PopNode_isClear(const PopNode *self);
void        PopNode_print(FILE * fp, PopNode * self, int indent);
void        PopNode_printShallow(PopNode * self, FILE * fp);
void       *PopNode_root(void * vself);
void        PopNode_sanityFromLeaf(PopNode * self, const char *file,
                                   int line);
int         PopNode_nsamples(PopNode * self);
void        PopNode_shiftParamPtrs(PopNode *self, size_t dp, int sign);
void        PopNode_shiftPopNodePtrs(PopNode *self, size_t dp, int sign);

NodeStore  *NodeStore_new(unsigned len, size_t elsize, void * v);
void        NodeStore_free(NodeStore *self);
void       *NodeStore_alloc(NodeStore * self);

#endif
