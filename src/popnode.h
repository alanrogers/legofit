#ifndef ARR_POPNODE_H
#  define ARR_POPNODE_H

#  include "typedefs.h"
#  include <stdbool.h>
#  include <gsl/gsl_rng.h>

int         PopNode_addChild(void * vparent, void * vchild);
void        PopNode_clear(PopNode * self);
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

#endif
