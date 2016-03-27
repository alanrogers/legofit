#ifndef ARR_GPTREE_H
#  define ARR_GPTREE_H

#  include "typedefs.h"
#  include <gsl/gsl_rng.h>

#  define POPNAMESIZE 30
#  define MAXSAMP ((int)(8*sizeof(tipId_t)))

void        Gene_tabulate(Gene * self, BranchTab * bt);
void        Gene_free(Gene * gene);
PopNode    *PopNode_new(double *twoNptr, double *tPtr);
void        PopNode_addChild(PopNode * parent, PopNode * child);
void        PopNode_mix(PopNode * child, double m, PopNode * introgressor,
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

#endif
