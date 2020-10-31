#ifndef ARR_POPNODE_H
#  define ARR_POPNODE_H

#  include "typedefs.h"
#  include <stdbool.h>
#  include <gsl/gsl_rng.h>

<<<<<<< HEAD
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

PopNode    *PopNode_new(double *twoN, double *start, NodeStore *ns);
int         PopNode_addChild(PopNode * parent, PopNode * child);
int         PopNode_mix(PopNode * child, double *mPtr, PopNode * introgressor,
                        PopNode * native);
void        PopNode_newGene(PopNode * self, unsigned ndx);
void        PopNode_addSample(PopNode * self, Gene * gene);
=======
int         PopNode_addChild(void * vparent, void * vchild);
void        PopNode_clear(PopNode * self);
>>>>>>> matcoal
Gene       *PopNode_coalesce(PopNode * self, gsl_rng * rng);
PopNode    *PopNode_dup(PopNode *old_root, PtrPtrMap *ppm);
int         PopNode_equals(PopNode *a, PopNode *b);
int         PopNode_feasible(const PopNode *self, Bounds bnd, int verbose);
void        PopNode_free(PopNode *self);
int         PopNode_isClear(const PopNode *self);
int         PopNode_mix(void * vchild, int mix_i, void * vintrogressor,
                        void * vnative, ParStore *ps);
void       *PopNode_new(int twoN_i, int start_i, ParStore *ps);
void        PopNode_newSample(PopNode * self, unsigned ndx);
void        PopNode_print(FILE * fp, void * vself, int indent);
void       *PopNode_root(void * vself);
void        PopNode_update(PopNode *self, ParStore *ps);
void        PopNode_unvisit(PopNode *self);
void        PopNode_shiftPopNodePtrs(PopNode *self, size_t dp, int sign);

<<<<<<< HEAD
NodeStore  *NodeStore_new(int len, PopNode *v);
void        NodeStore_free(NodeStore *self);
PopNode    *NodeStore_alloc(NodeStore *self);

=======
>>>>>>> matcoal
#endif
