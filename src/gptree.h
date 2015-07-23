#ifndef ARR_GPTREE_H
#define ARR_GPTREE_H

#include "typedefs.h"
#include <gsl/gsl_rng.h>

#define POPNAMESIZE 30
#define MAXSAMP ((int)(8*sizeof(tipId_t)))

Gene *Gene_new(tipId_t tipId);
void  Gene_tabulate(Gene *self, BranchTab *bt);
void  Gene_addToBranch(Gene *gene, double x);
Gene *Gene_join(Gene *lchild, Gene *rchild);
void  Gene_free(Gene *gene);
double Gene_lastInterval(Gene *gene, tipId_t *tipId);
double Gene_checkInterval(Gene *gene, tipId_t *tipId, double *branch);
double Gene_getRightLen(Gene *gene, tipId_t tipId);
PopNode *PopNode_new(double twoN, double start);
void PopNode_set(PopNode *self, double twoN, double start, double end);
void PopNode_addChild(PopNode *parent, PopNode *child);
void PopNode_mix(PopNode *child, double m, PopNode *introgressor, PopNode *native);
void PopNode_join(PopNode *parent, PopNode *lchild, PopNode *rchild);
void PopNode_newGene(PopNode *pnode, unsigned ndx);
void PopNode_addSample(PopNode *pnode, Gene *gene);
Gene * PopNode_coalesce(PopNode *pnode, gsl_rng *rng);
void PopNode_free(PopNode *pnode);
void PopNode_clear(PopNode *pnode);
void PopNode_print(FILE *fp, PopNode *pnode, int indent);
void PopNode_printShallow(PopNode *self, FILE *fp);
PopNode *PopNode_root(PopNode *self);
void PopNode_sanityFromLeaf(PopNode *self, const char *file, int line);
int PopNode_nsamples(PopNode *self);
double survival(double t, double twoN);

#endif
