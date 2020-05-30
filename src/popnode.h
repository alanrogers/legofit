#ifndef ARR_POPNODE_H
#  define ARR_POPNODE_H

#  include "typedefs.h"
#  include <stdbool.h>
#  include <gsl/gsl_rng.h>

struct PopNode {
    int         nparents, nchildren, nsamples;
    double      *twoN;           // ptr to current pop size
    double      *start, *end;    // duration of this PopNode
    double      *mix;            // ptr to frac of pop derived from parent[1]
    struct PopNode *parent[2];
    struct PopNode *child[2];

    Gene       *sample[MAXSAMP]; // not locally owned
};

int         PopNode_addChild(void * vparent, void * vchild);
void        PopNode_clear(PopNode * self);
Gene       *PopNode_coalesce(PopNode * self, gsl_rng * rng);
int         PopNode_feasible(const PopNode *self, Bounds bnd, int verbose);
int         PopNode_isClear(const PopNode *self);
int         PopNode_mix(void * vchild, double *mPtr, void * vintrogressor,
                        void * vnative);
void       *PopNode_new(double *twoN, double *start, NodeStore *ns);
void        PopNode_newGene(PopNode * self, unsigned ndx);
void        PopNode_print(FILE * fp, void * vself, int indent);
void       *PopNode_root(void * vself);
void        PopNode_shiftParamPtrs(PopNode *self, size_t dp, int sign);
void        PopNode_shiftPopNodePtrs(PopNode *self, size_t dp, int sign);

#endif
