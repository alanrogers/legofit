#ifndef ARR_GENE
#define ARR_GENE

#include "branchtab.h"

struct Gene {
    tipId_t     tipId;
    struct Gene *parent, *lchild, *rchild;
    double      branch;
};

Gene       *Gene_join(Gene * lchild, Gene * rchild);
Gene       *Gene_new(tipId_t tipId);
void        Gene_tabulate(Gene * self, BranchTab * bt);
void        Gene_free(Gene * gene);

static inline void Gene_addToBranch(Gene * gene, double x);

static inline void Gene_addToBranch(Gene * gene, double x) {
    gene->branch += x;
}


#endif
