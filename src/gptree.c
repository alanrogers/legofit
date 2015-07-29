/**
 * @file gptree.c
 * @brief Methods for simulating gene genealogies within a given tree
 * of populations, and allowing populations to mix and also to split.
 */

#include "gptree.h"
#include "misc.h"
#include "binary.h"
#include "branchtab.h"
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

struct Gene {
    tipId_t     tipId;
    struct Gene *parent, *lchild, *rchild;
    double      branch;
};

struct PopNode {
    int         nparents, nchildren, nsamples;
    double      twoN;           // current pop size to ancestral pop size 
    double      start, end;     // duration of this PopNode
    double      mix;            // mix=frac of pop derived from parent[1]
    struct PopNode *parent[2];
    struct PopNode *child[2];
    Gene       *sample[MAXSAMP];
};

static void PopNode_sanityCheck(PopNode * self, const char *file, int lineno);
static inline void Gene_addToBranch(Gene * gene, double x);
static Gene *Gene_join(Gene * lchild, Gene * rchild);
static Gene *Gene_new(tipId_t tipId);

static Gene *Gene_new(tipId_t tipId) {
    Gene       *gene = malloc(sizeof(Gene));
    checkmem(gene, __FILE__, __LINE__);

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
void Gene_tabulate(Gene * self, BranchTab * bt) {
    if(self == NULL)
        return;

    if(self->parent && !isPow2(self->tipId))
        BranchTab_add(bt, self->tipId, self->branch);

    Gene_tabulate(self->lchild, bt);
    Gene_tabulate(self->rchild, bt);
}

static Gene *Gene_join(Gene * lchild, Gene * rchild) {
    tipId_t     id = lchild->tipId | rchild->tipId;
    Gene       *parent = Gene_new(id);
    checkmem(parent, __FILE__, __LINE__);
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

static inline void Gene_addToBranch(Gene * gene, double x) {
    gene->branch += x;
}

void PopNode_sanityFromLeaf(PopNode * self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self != NULL, file, line);
    switch (self->nparents) {
    case 0:
        REQUIRE(self->parent[0] == NULL, file, line);
        REQUIRE(self->parent[1] == NULL, file, line);
        REQUIRE(self->mix == 0.0, file, line);
        if(!isinf(self->end)) {
            fflush(stdout);
            PopNode_printShallow(self, stderr);
        }
        REQUIRE(isinf(self->end), file, line);
        break;
    case 1:
        REQUIRE(self->parent[0] != NULL, file, line);
        REQUIRE(self->parent[1] == NULL, file, line);
        if(self->mix != 0.0)
            PopNode_printShallow(self, stdout);
        REQUIRE(self->mix == 0.0, file, line);
        break;
    default:
        REQUIRE(self->nparents == 2, file, line);
        REQUIRE(self->parent[0] != NULL, file, line);
        REQUIRE(self->parent[1] != NULL, file, line);
        REQUIRE(self->mix >= 0.0, file, line);
        break;
    }
    switch (self->nchildren) {
    case 0:
        REQUIRE(self->child[0] == NULL, file, line);
        REQUIRE(self->child[1] == NULL, file, line);
        break;
    case 1:
        REQUIRE(self->child[0] != NULL, file, line);
        REQUIRE(self->child[1] == NULL, file, line);
        break;
    default:
        REQUIRE(self->nchildren == 2, file, line);
        REQUIRE(self->child[0] != NULL, file, line);
        REQUIRE(self->child[1] != NULL, file, line);
        break;
    }
    REQUIRE(self->start < self->end, file, line);
    if(self->nparents > 0)
        PopNode_sanityFromLeaf(self->parent[0], file, line);
    if(self->nparents > 1)
        PopNode_sanityFromLeaf(self->parent[1], file, line);
#endif
}

/// Find root of population tree, starting from given node.
PopNode    *PopNode_root(PopNode * self) {
    PopNode    *r0, *r1;
    assert(self);
    switch (self->nparents) {
    case 0:
        return self;
        break;
    case 1:
        return PopNode_root(self->parent[0]);
        break;
    case 2:
        r0 = PopNode_root(self->parent[0]);
        r1 = PopNode_root(self->parent[1]);
        if(r0 != r1) {
            fprintf(stderr, "%s:%s:%d: Node has multiple roots\n",
                    __FILE__, __func__, __LINE__);
            exit(EXIT_FAILURE);
        }
        return r0;
        break;
    default:
        fprintf(stderr, "%s:%s:%d: Node %d parents\n",
                __FILE__, __func__, __LINE__, self->nparents);
        exit(EXIT_FAILURE);
    }
    /* NOTREACHED */
    return NULL;
}

/** Remove all references to samples from tree of populations */
void PopNode_clear(PopNode * self) {
    int         i;
    for(i = 0; i < self->nchildren; ++i)
        PopNode_clear(self->child[i]);

    self->nsamples = 0;
    memset(self->sample, 0, sizeof(self->sample));
    PopNode_sanityCheck(self, __FILE__, __LINE__);
}

void PopNode_print(FILE * fp, PopNode * self, int indent) {
    int         i;
    for(i = 0; i < indent; ++i)
        fputs("   ", fp);
    fprintf(fp, "%p twoN=%lf ntrval=(%lf,", self, self->twoN, self->start);
    if(self->end < DBL_MAX)
        fprintf(fp, "%lf)\n", self->end);
    else
        fprintf(fp, "Inf)\n");

    for(i = 0; i < self->nchildren; ++i)
        PopNode_print(fp, self->child[i], indent + 1);
}

void PopNode_printShallow(PopNode * self, FILE * fp) {
    fprintf(fp, "%p twoN=%lf mix=%lf ntrval=(%lf,",
            self, self->twoN, self->mix, self->start);
    if(self->end < DBL_MAX)
        fprintf(fp, "%lf)", self->end);
    else
        fprintf(fp, "Inf)");

    switch (self->nparents) {
    case 0:
        fprintf(fp, " par=0");
        break;
    case 1:
        fprintf(fp, " par=%p", self->parent[0]);
        break;
    default:
        fprintf(fp, " par=[%p,%p]", self->parent[0], self->parent[1]);
        break;
    }
    switch (self->nchildren) {
    case 0:
        fprintf(fp, " child=0");
        break;
    case 1:
        fprintf(fp, " child=%p", self->child[0]);
        break;
    default:
        fprintf(fp, " child=[%p,%p]", self->child[0], self->child[1]);
        break;
    }
    putc('\n', fp);
}

int PopNode_nsamples(PopNode * self) {
    return self->nsamples;
}

PopNode    *PopNode_new(double twoN, double start) {
    PopNode    *new = malloc(sizeof(PopNode));
    checkmem(new, __FILE__, __LINE__);

    new->nparents = new->nchildren = new->nsamples = 0;
    new->twoN = twoN;
    new->mix = 0.0;
    new->start = start;
    new->end = HUGE_VAL;

    memset(new->sample, 0, sizeof(new->sample));
    memset(new->parent, 0, sizeof(new->parent));
    memset(new->child, 0, sizeof(new->child));

    PopNode_sanityCheck(new, __FILE__, __LINE__);
    return new;
}

/// Connect parent and child 
void PopNode_addChild(PopNode * parent, PopNode * child) {
    if(parent->nchildren > 1)
        eprintf("%s:%s:%d: Can't add child because parent already has %d.\n",
                __FILE__, __func__, __LINE__, parent->nchildren);
    if(child->nparents > 1)
        eprintf("%s:%s:%d: Can't add parent because child already has %d.\n",
                __FILE__, __func__, __LINE__, child->nparents);
    if(child->start > parent->start)
        eprintf("%s:%s:%d: Child start (%lf) must be < parent start (%lf)\n",
                __FILE__, __func__, __LINE__, child->start, parent->start);
    if(child->end < HUGE_VAL) {
        if(child->end != parent->start)
            eprintf
                ("%s:%s:%d: Date mismatch. child->end=%lf != %lf = parent->start\n",
                 __FILE__, __func__, __LINE__, child->end, parent->start);
    } else
        child->end = parent->start;
    parent->child[parent->nchildren] = child;
    child->parent[child->nparents] = parent;
    ++parent->nchildren;
    ++child->nparents;
    PopNode_sanityCheck(parent, __FILE__, __LINE__);
    PopNode_sanityCheck(child, __FILE__, __LINE__);
}

static void PopNode_sanityCheck(PopNode * self, const char *file, int lineno) {
#ifndef NDEBUG
    int         i;

    REQUIRE(self != NULL, file, lineno);

    for(i = 0; i < self->nsamples; ++i)
        REQUIRE(self->sample[i] != NULL, file, lineno);
#endif
}

void PopNode_addSample(PopNode * self, Gene * gene) {
    if(self->nsamples == MAXSAMP) {
        fprintf(stderr, "%s:%s:%d: Too many samples\n",
                __FILE__, __func__, __LINE__);
        exit(1);
    }
    self->sample[self->nsamples] = gene;
    ++self->nsamples;
    PopNode_sanityCheck(self, __FILE__, __LINE__);
}

void PopNode_mix(PopNode * child, double m, PopNode * introgressor,
                 PopNode * native) {
    if(introgressor->nchildren > 1)
        EPRINTF(("%s:%s:%d: Can't add child because introgressor already has %d.\n", __FILE__, __func__, __LINE__, introgressor->nchildren));
    if(native->nchildren > 1)
        EPRINTF(("%s:%s:%d: Can't add child because native parent already has %d.\n", __FILE__, __func__, __LINE__, native->nchildren));
    if(child->nparents > 0)
        EPRINTF(("%s:%s:%d: Can't add 2 parents because child already has %d.\n", __FILE__, __func__, __LINE__, child->nparents));
    if(child->end < HUGE_VAL) {
        if(child->end != introgressor->start)
            EPRINTF(("%s:%s:%d: Date mismatch. child->end=%lf != %lf = introgressor->start\n", __FILE__, __func__, __LINE__, child->end, introgressor->start));
        if(child->end != native->start)
            EPRINTF(("%s:%s:%d: Date mismatch. child->end=%lf != %lf = native->start\n", __FILE__, __func__, __LINE__, child->end, native->start));
    } else
        child->end = native->start;

    child->parent[0] = native;
    child->parent[1] = introgressor;
    child->nparents = 2;
    child->mix = m;
    introgressor->child[introgressor->nchildren] = child;
    ++introgressor->nchildren;
    native->child[native->nchildren] = child;
    ++native->nchildren;
    PopNode_sanityCheck(child, __FILE__, __LINE__);
    PopNode_sanityCheck(introgressor, __FILE__, __LINE__);
    PopNode_sanityCheck(native, __FILE__, __LINE__);
}

void PopNode_newGene(PopNode * self, unsigned ndx) {
    assert(1 + self->nsamples < MAXSAMP);

    assert(ndx < sizeof(tipId_t));
    Gene       *gene = Gene_new(1UL << ndx);
    checkmem(gene, __FILE__, __LINE__);
    self->sample[self->nsamples] = gene;
    ++self->nsamples;
    PopNode_sanityCheck(self, __FILE__, __LINE__);
}

/// Coalesce gene tree within population tree.
Gene       *PopNode_coalesce(PopNode * self, gsl_rng * rng) {
    unsigned long i, j, k;
    double      x;

    if(self->child[0])
        (void) PopNode_coalesce(self->child[0], rng);
    if(self->child[1])
        (void) PopNode_coalesce(self->child[1], rng);

    double      t = self->start;
    if(!(t < self->end))
        PopNode_print(stdout, self, 0);
    assert(t < self->end);

    // Coalescent loop continues until only one sample is left
    // or we reach the end of the interval.
    while(self->nsamples > 1 && t < self->end) {
        double      m;
        {
            int         n = self->nsamples;
            m = 2.0 * self->twoN / (n * (n - 1));
        }
        x = gsl_ran_exponential(rng, m);

        if(t + x < self->end) {
            // coalescent event within interval
            t += x;
            for(i = 0; i < self->nsamples; ++i)
                Gene_addToBranch(self->sample[i], x);

            // choose a random pair to join
            i = gsl_rng_uniform_int(rng, self->nsamples);
            j = gsl_rng_uniform_int(rng, self->nsamples - 1);
            if(j >= i)
                ++j;
            if(j < i) {
                k = i;
                i = j;
                j = k;
            }
            assert(i < j);

            self->sample[i] = Gene_join(self->sample[i], self->sample[j]);
            checkmem(self->sample[i], __FILE__, __LINE__);
            --self->nsamples;
            if(j != self->nsamples) {
                self->sample[j] = self->sample[self->nsamples];
                self->sample[self->nsamples] = NULL;
            }
        } else {
            // no coalescent event within interval
            x = self->end - t;
            for(i = 0; i < self->nsamples; ++i)
                Gene_addToBranch(self->sample[i], x);
            t = self->end;
        }
    }

    // Make sure we're at the end of the interval
    if(t < self->end) {
        assert(self->nsamples < 2);
        x = self->end - t;
        for(i = 0; i < self->nsamples; ++i)
            Gene_addToBranch(self->sample[i], x);
        t = self->end;
    }
    // If we have both samples and parents, then move samples to parents
    if(self->nsamples > 0 && self->nparents > 0) {
        assert(t == self->end);
        for(i = 0; i < self->nsamples; ++i) {
            // move current sample to parental population 
            if(self->nparents > 1 && gsl_rng_uniform(rng) < self->mix) {
                PopNode_addSample(self->parent[1], self->sample[i]);
            } else if(self->nparents > 0)
                PopNode_addSample(self->parent[0], self->sample[i]);
        }
        self->nsamples = 0;
    }

    PopNode_sanityCheck(self, __FILE__, __LINE__);
    return (self->nsamples == 1 ? self->sample[0] : NULL);
}

/// Free node but not descendants
void PopNode_free(PopNode * self) {
    free(self);
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
            fprintf(stderr,"usage: xgptree [-v]\n");
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
    Gene_tabulate(g5, bt);

    assert(BranchTab_size(bt) == 1);
    assert(2.0 == BranchTab_get(bt, (id1|id2)));

    unitTstResult("Gene", "OK");

    double twoN0 = 1.0, t0= 0.0;
    PopNode *p0 = PopNode_new(twoN0, t0);
    assert(p0->twoN == twoN0);
    assert(p0->start == t0);
    assert(p0->mix == 0.0);
    assert(p0->end == HUGE_VAL);
    assert(p0->nsamples == 0);
    assert(p0->nchildren == 0);
    assert(p0->child[0] == NULL);
    assert(p0->child[1] == NULL);
    assert(p0->parent[0] == NULL);
    assert(p0->parent[1] == NULL);

    double twoN1 = 100.0, t1= 123.0;
    PopNode *p1 = PopNode_new(twoN1, t1);
    assert(p1->twoN == twoN1);
    assert(p1->start == t1);
    assert(p1->mix == 0.0);
    assert(p1->end == HUGE_VAL);
    assert(p1->nsamples == 0);
    assert(p1->nchildren == 0);
    assert(p1->child[0] == NULL);
    assert(p1->child[1] == NULL);
    assert(p1->parent[0] == NULL);
    assert(p1->parent[1] == NULL);

    g1 = Gene_new(id1);
    g2 = Gene_new(id2);
    PopNode_addSample(p1, g1);
    PopNode_addSample(p1, g2);
    assert(p1->nsamples == 2);

    PopNode_addChild(p1, p0);
    assert(p1->nchildren == 1);
    assert(p0->nparents == 1);
    assert(p1->child[0] == p0);
    assert(p0->parent[0] == p1);

    unitTstResult("PopNode", "untested");

    return 0;
}
#endif
