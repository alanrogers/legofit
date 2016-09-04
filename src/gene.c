/**
 * @file gene.c
 * @brief Class Gene. Defines the objects that are linked together to
 * form a gene genealogy.
 */

#include "gene.h"
#include "misc.h"
#include "binary.h"
#include <string.h>
#include <stdlib.h>

Gene *Gene_new(tipId_t tipId) {
    Gene       *gene = malloc(sizeof(Gene));
    CHECKMEM(gene);

    gene->tipId = tipId;
    gene->parent = gene->lchild = gene->rchild = NULL;
    gene->branch = 0.0;

    return gene;
}

/**
 * Tabulate branch lengths by tipId.
 *
 * Ignore root (which has null parent) because mutations
 * there don't contribute to genetic variation.
 * 
 * Ignore nodes with only one descendant, which are recognizable because 
 * their tipIds are powers of 2.
 * 
 * Tabulate everything else.
 */
void Gene_tabulate(Gene * self, BranchTab * bt, int doSing) {
    if(self == NULL)
        return;

    if(self->parent && (doSing || !isPow2(self->tipId)))
        BranchTab_add(bt, self->tipId, self->branch);

    Gene_tabulate(self->lchild, bt, doSing);
    Gene_tabulate(self->rchild, bt, doSing);
}

Gene *Gene_join(Gene * lchild, Gene * rchild) {
    tipId_t     id = lchild->tipId | rchild->tipId;
    Gene       *parent = Gene_new(id);
    CHECKMEM(parent);
    parent->lchild = lchild;
    parent->rchild = rchild;
    lchild->parent = rchild->parent = parent;
    return parent;
}

void Gene_free(Gene * gene) {
    if(gene == NULL)
        return;
    Gene_free(gene->lchild);
    Gene_free(gene->rchild);
    free(gene);
}

#ifdef TEST

#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {

    int verbose=0;

    if(argc > 1) {
        if(argc!=2 || 0!=strcmp(argv[1], "-v")) {
            fprintf(stderr,"usage: xgene [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    tipId_t id1 = 1;
    Gene *g1 = Gene_new(id1);
    assert(g1);
    Gene_addToBranch(g1, 1.0);
    assert(g1->tipId == id1);
    assert(g1->parent == NULL);
    assert(g1->lchild == NULL);
    assert(g1->rchild == NULL);
    assert(g1->branch == 1.0);

    tipId_t id2 = 2;
    Gene *g2 = Gene_new(id2);
    assert(g2);
    Gene_addToBranch(g2, 1.0);
    assert(g2->tipId == id2);
    assert(g2->parent == NULL);
    assert(g2->lchild == NULL);
    assert(g2->rchild == NULL);
    assert(g2->branch == 1.0);

    Gene *g3 = Gene_join(g1, g2);
    assert(g3);
    Gene_addToBranch(g3, 2.0);
    assert(g3->tipId == (id1|id2));
    assert(g3->parent == NULL);
    assert(g3->lchild == g1);
    assert(g3->rchild == g2);
    assert(g3->branch == 2.0);

    tipId_t id4 = 4;
    Gene *g4 = Gene_new(id4);
    assert(g4);
    Gene_addToBranch(g4, 2.0);
    assert(g4->tipId == id4);
    assert(g4->parent == NULL);
    assert(g4->lchild == NULL);
    assert(g4->rchild == NULL);
    assert(g4->branch == 2.0);

    Gene *g5 = Gene_join(g3, g4);
    assert(g5);
    assert(g5->tipId == (id1|id2|id4));
    assert(g5->parent == NULL);
    assert(g5->lchild == g3);
    assert(g5->rchild == g4);
    assert(g5->branch == 0.0);

    BranchTab *bt = BranchTab_new();
    Gene_tabulate(g5, bt, 0);

    assert(BranchTab_size(bt) == 1);
    assert(2.0 == BranchTab_get(bt, (id1|id2)));

    unitTstResult("Gene", "OK");

    return 0;
}
#endif
