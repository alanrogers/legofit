#ifndef ARR_POPNODE_H
#  define ARR_POPNODE_H

#  include "typedefs.h"
#  include <stdbool.h>
#  include <gsl/gsl_rng.h>

struct PopNode {
    int         visited; // has the coalescent visited this node yet?
    int         nparents, nchildren, nsamples;
    double      twoN;            // haploid pop size
    double      start, end;      // duration of this PopNode
    double      mix;             // frac of pop derived from parent[1]

    // indices into ParStore array
    int twoN_i, start_i, end_i, mix_i;

    struct PopNode *parent[2];
    struct PopNode *child[2];

    Gene       *sample[MAXSAMP]; // not locally owned
};

int         PopNode_addChild(void * vparent, void * vchild);
void        PopNode_clear(PopNode * self);
Gene       *PopNode_coalesce(PopNode * self, gsl_rng * rng);
int         PopNode_feasible(const PopNode *self, Bounds bnd, int verbose);
int         PopNode_isClear(const PopNode *self);
int         PopNode_mix(void * vchild, int mix_i, void * vintrogressor,
                        void * vnative, ParStore *ps);
void       *PopNode_new(int twoN_i, int start_i, ParStore *ps, NodeStore * ns);
void        PopNode_newGene(PopNode * self, unsigned ndx);
void        PopNode_print(FILE * fp, void * vself, int indent);
void       *PopNode_root(void * vself);
void        PopNode_update(PopNode *self, ParStore *ps);
void        PopNode_unvisit(PopNode *self);
void        PopNode_shiftPopNodePtrs(PopNode *self, size_t dp, int sign);

#endif
