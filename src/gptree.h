#ifndef ARR_GPTREE_H
#define ARR_GPTREE_H

typedef struct PopNode PopNode;
typedef struct Gene Gene;
typedef unsigned long tipId_t;

enum { MAXSAMPLES = 4 };

#include <gsl/gsl_rng.h>

Gene *Gene_new(tipId_t tipId);
void Gene_addToBranch(Gene *gene, double x);
Gene *Gene_join(Gene *lchild, Gene *rchild);
void Gene_free(Gene *gene);
double Gene_lastInterval(Gene *gene, tipId_t *tipId);
double Gene_checkInterval(Gene *gene, tipId_t *tipId, double *branch);
double Gene_getRightLen(Gene *gene, tipId_t tipId);
PopNode *PopNode_new(double K, double start, double end);
void PopNode_set(PopNode *self, double K, double start, double end);
void PopNode_addChild(PopNode *parent, PopNode *child);
void PopNode_mix(PopNode *pnode, double m, PopNode *immigrant, PopNode *native);
void PopNode_endToEnd(PopNode *pnode, PopNode *ancestor);
void PopNode_join(PopNode *parent, PopNode *lchild, PopNode *rchild);
void PopNode_newGene(PopNode *pnode, unsigned ndx);
void PopNode_addSample(PopNode *pnode, Gene *gene);
Gene * PopNode_coalesce(PopNode *pnode, gsl_rng *rng);
void PopNode_free(PopNode *pnode);
void PopNode_clear(PopNode *pnode);
void PopNode_print(FILE *fp, PopNode *pnode, int indent);
PopNode *PopNode_root(PopNode *self);
double survival(double t, double K);

#endif
