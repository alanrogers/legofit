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
    tipId_t tipId;
    struct Gene *parent, *lchild, *rchild;
    double branch;
};

struct PopNode {
    int nparents, nchildren, nsamples;
    double twoN;       // current pop size to ancestral pop size 
    double start, end; // duration of this PopNode
    double mix;        // mix=frac of pop derived from parent[1]
    struct PopNode *parent[2];
    struct PopNode *child[2];
    Gene *sample[MAXSAMP];
};

static void PopNode_sanityCheck(PopNode *pnode, const char *file, int lineno);

void PopNode_sanityFromLeaf(PopNode *self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self!=NULL, file, line);
    switch(self->nparents) {
    case 0:
        REQUIRE(self->parent[0]==NULL, file, line);
        REQUIRE(self->parent[1]==NULL, file, line);
        REQUIRE(self->mix == 0.0, file, line);
        if(!isinf(self->end)) {
            fflush(stdout);
            PopNode_printShallow(self, stderr);
        }
        REQUIRE(isinf(self->end), file, line);
        break;
    case 1:
        REQUIRE(self->parent[0]!=NULL, file, line);
        REQUIRE(self->parent[1]==NULL, file, line);
        if(self->mix != 0.0)
            PopNode_printShallow(self, stdout);
        REQUIRE(self->mix == 0.0, file, line);
        break;
    default:
        REQUIRE(self->nparents==2, file, line);
        REQUIRE(self->parent[0]!=NULL, file, line);
        REQUIRE(self->parent[1]!=NULL, file, line);
        REQUIRE(self->mix >= 0.0, file, line);
        break;
    }
    switch(self->nchildren) {
    case 0:
        REQUIRE(self->child[0]==NULL, file, line);
        REQUIRE(self->child[1]==NULL, file, line);
        break;
    case 1:
        REQUIRE(self->child[0]!=NULL, file, line);
        REQUIRE(self->child[1]==NULL, file, line);
        break;
    default:
        REQUIRE(self->nchildren==2, file, line);
        REQUIRE(self->child[0]!=NULL, file, line);
        REQUIRE(self->child[1]!=NULL, file, line);
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
PopNode *PopNode_root(PopNode *self) {
    PopNode *r0, *r1;
    assert(self);
    switch(self->nparents) {
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
            fprintf(stderr,"%s:%s:%d: Node has multiple roots\n",
                    __FILE__, __func__, __LINE__);
            exit(EXIT_FAILURE);
        }
        return r0;
        break;
    default:
            fprintf(stderr,"%s:%s:%d: Node %d parents\n",
                    __FILE__, __func__, __LINE__, self->nparents);
            exit(EXIT_FAILURE);
    }
    /* NOTREACHED */
    return NULL;
}

void PopNode_set(PopNode *self, double twoN, double start, double end) {
    self->twoN = twoN;
    self->start = start;
    self->end = end;
}

/** Remove all references to samples from tree of populations */
void PopNode_clear(PopNode *pnode) {
    int i;
    for(i=0; i < pnode->nchildren; ++i)
        PopNode_clear(pnode->child[i]);

    pnode->nsamples = 0;
    memset(pnode->sample, 0, sizeof(pnode->sample));
    PopNode_sanityCheck(pnode, __FILE__, __LINE__);
}

void PopNode_print(FILE *fp, PopNode *pnode, int indent) {
    int i;
    for(i=0; i< indent; ++i)
        fputs("   ", fp);
    fprintf(fp, "%p twoN=%lf ntrval=(%lf,",
            pnode, pnode->twoN, pnode->start);
    if(pnode->end < DBL_MAX)
        fprintf(fp, "%lf)\n", pnode->end);
    else
        fprintf(fp, "Inf)\n");

    for(i=0; i < pnode->nchildren; ++i)
        PopNode_print(fp, pnode->child[i], indent+1);
}

void PopNode_printShallow(PopNode *self, FILE *fp) {
    fprintf(fp, "%p twoN=%lf mix=%lf ntrval=(%lf,",
            self, self->twoN, self->mix, self->start);
    if(self->end < DBL_MAX)
        fprintf(fp, "%lf)", self->end);
    else
        fprintf(fp, "Inf)");

    switch(self->nparents) {
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
    switch(self->nchildren) {
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

int PopNode_nsamples(PopNode *self) {
    return self->nsamples;
}

Gene *Gene_new(tipId_t tipId) {
    Gene *gene = malloc(sizeof(Gene));
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
void Gene_tabulate(Gene *self, BranchTab *bt) {
    if(self==NULL)
        return;

    if(self->parent && !isPow2(self->tipId))
        BranchTab_add(bt, self->tipId, self->branch);

    Gene_tabulate(self->lchild, bt);
    Gene_tabulate(self->rchild, bt);
}

void Gene_addToBranch(Gene *gene, double x) {
    gene->branch += x;
}

Gene *Gene_join(Gene *lchild, Gene *rchild) {
    tipId_t id = lchild->tipId | rchild->tipId;
    Gene *parent = Gene_new(id);
    checkmem(parent, __FILE__, __LINE__);
    parent->lchild = lchild;
    parent->rchild = rchild;
    lchild->parent = rchild->parent = parent;
    return parent;
}

/**
 * Return the length of the deepest coalescent interval.
 * Set *tipId equal to the tipId of the oldest of the two children
 * of the given node.
 */
double Gene_lastInterval(Gene *gene, tipId_t *tipId) {
    if(gene->lchild==NULL || gene->rchild==NULL)
        eprintf("%s:%s%d: root has null children:"
                " can't calculate branch length\n",
                __FILE__,__func__,__LINE__);
    Gene *child;
    if(gene->lchild->branch < gene->rchild->branch)
        child = gene->lchild;
    else
        child = gene->rchild;
    *tipId = child->tipId;
    return child->branch;
}

double Gene_checkInterval(Gene *gene, tipId_t *tipId, double *branch) {
    if(gene->lchild==NULL || gene->rchild==NULL)
	return 0.0;
    /**
     * Check if either of the node's children are the target. The final
     * interval can never be the target, so it is safe to assume that if
     * we are at this point we don't need to check equality on the current
     * gene.
     */
    if(gene->lchild->tipId == *tipId)
	    *branch = gene->lchild->branch;
    else if(gene->rchild->tipId == *tipId)
	    *branch = gene->rchild->branch;
    else {
	    Gene_checkInterval(gene->lchild, tipId, branch);
	    Gene_checkInterval(gene->rchild, tipId, branch);
    }
    if (branch == NULL)
	    return 0.0;
    else
	    return *branch;
}

double Gene_getRightLen(Gene *gene, tipId_t tipId) {
	if(gene == NULL)
		return 0.0;
	if(gene->tipId == tipId)
		return gene->branch;
    return fmax(Gene_getRightLen(gene->lchild, tipId),
				Gene_getRightLen(gene->rchild, tipId));
}

void Gene_free(Gene *gene) {
    if(gene == NULL)
        return;
    Gene_free(gene->lchild);
    Gene_free(gene->rchild);
    free(gene);
}

PopNode *PopNode_new(double twoN, double start) {
    PopNode *pnode = malloc(sizeof(PopNode));
    checkmem(pnode, __FILE__, __LINE__);

    pnode->nparents = pnode->nchildren = pnode->nsamples = 0;
    pnode->twoN = twoN;
    pnode->mix = 0.0;
    pnode->start = start;
    pnode->end = HUGE_VAL;

    memset(pnode->sample, 0, sizeof(pnode->sample));
    memset(pnode->parent, 0, sizeof(pnode->parent));
    memset(pnode->child, 0, sizeof(pnode->child));

    PopNode_sanityCheck(pnode, __FILE__, __LINE__);
    return pnode;
}

/// Connect parent and child 
void PopNode_addChild(PopNode *parent, PopNode *child) {
    if(parent->nchildren > 1)
        eprintf("%s:%s:%d: Can't add child because parent already has %d.\n",
                 __FILE__,__func__,__LINE__, parent->nchildren);
    if(child->nparents > 1)
        eprintf("%s:%s:%d: Can't add parent because child already has %d.\n",
                 __FILE__,__func__,__LINE__, child->nparents);
    if(child->start > parent->start)
        eprintf("%s:%s:%d: Child start (%lf) must be < parent start (%lf)\n",
                __FILE__, __func__, __LINE__,
                child->start, parent->start);
    if(child->end < HUGE_VAL) {
        if(child->end != parent->start)
            eprintf("%s:%s:%d: Date mismatch. child->end=%lf != %lf = parent->start\n",
                     __FILE__,__func__,__LINE__,
                     child->end, parent->start);
    }else
        child->end = parent->start;
    parent->child[parent->nchildren] = child;
    child->parent[child->nparents] = parent;
    ++parent->nchildren;
    ++child->nparents;
    PopNode_sanityCheck(parent, __FILE__, __LINE__);
    PopNode_sanityCheck(child, __FILE__, __LINE__);
}

static void PopNode_sanityCheck(PopNode *pnode, const char *file, int lineno) {
#ifndef NDEBUG
    int i;

    REQUIRE(pnode != NULL, file, lineno);

    for(i=0; i < pnode->nsamples; ++i)
        REQUIRE(pnode->sample[i] != NULL, file, lineno);
#endif
}

void PopNode_addSample(PopNode *pnode, Gene *gene) {
    if(pnode->nsamples == MAXSAMP) {
        fprintf(stderr,"%s:%s:%d: Too many samples\n",
                __FILE__,__func__,__LINE__);
        exit(1);
    }
    pnode->sample[pnode->nsamples] = gene;
    ++pnode->nsamples;
    PopNode_sanityCheck(pnode, __FILE__, __LINE__);
}

void PopNode_mix(PopNode *child, double m, PopNode *introgressor, PopNode *native) {
    if(introgressor->nchildren > 1)
        EPRINTF(("%s:%s:%d: Can't add child because introgressor already has %d.\n",
                 __FILE__,__func__,__LINE__, introgressor->nchildren));
    if(native->nchildren > 1)
        EPRINTF(("%s:%s:%d: Can't add child because native parent already has %d.\n",
                 __FILE__,__func__,__LINE__, native->nchildren));
    if(child->nparents > 0)
        EPRINTF(("%s:%s:%d: Can't add 2 parents because child already has %d.\n",
                 __FILE__,__func__,__LINE__, child->nparents));
    if(child->end < HUGE_VAL) {
        if(child->end != introgressor->start)
            EPRINTF(("%s:%s:%d: Date mismatch. child->end=%lf != %lf = introgressor->start\n",
                     __FILE__,__func__,__LINE__,
                     child->end, introgressor->start));
        if(child->end != native->start)
            EPRINTF(("%s:%s:%d: Date mismatch. child->end=%lf != %lf = native->start\n",
                     __FILE__,__func__,__LINE__,
                     child->end, native->start));
    }else
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

void PopNode_join(PopNode *parent, PopNode *lchild, PopNode *rchild) {
    PopNode_addChild(parent, lchild);
    PopNode_addChild(parent, rchild);
    lchild->nparents = rchild->nparents = 1;
    lchild->parent[0] = rchild->parent[0] = parent;
    lchild->mix = rchild->mix = 0.0;
    PopNode_sanityCheck(parent, __FILE__, __LINE__);
    PopNode_sanityCheck(lchild, __FILE__, __LINE__);
    PopNode_sanityCheck(rchild, __FILE__, __LINE__);
}

void PopNode_newGene(PopNode *pnode, unsigned ndx) {
    assert(1 + pnode->nsamples < MAXSAMP);

    assert(ndx < sizeof(tipId_t));
    Gene *gene = Gene_new(1UL << ndx);
    checkmem(gene, __FILE__, __LINE__);
    pnode->sample[pnode->nsamples] = gene;
    ++pnode->nsamples;
    PopNode_sanityCheck(pnode, __FILE__, __LINE__);
}

/// Coalesce gene tree within population tree.
Gene * PopNode_coalesce(PopNode *pnode, gsl_rng *rng) {
    unsigned long i, j, k;
    double x;

    if(pnode->child[0])
        (void) PopNode_coalesce(pnode->child[0], rng);
    if(pnode->child[1])
        (void) PopNode_coalesce(pnode->child[1], rng);

    double t = pnode->start;
    if(!(t < pnode->end))
        PopNode_print(stdout, pnode, 0);
    assert(t < pnode->end);

    // Coalescent loop continues until only one sample is left
    // or we reach the end of the interval.
    while(pnode->nsamples > 1 && t < pnode->end) {
        double m;
        {
            int n = pnode->nsamples;
            m = 2.0 * pnode->twoN / (n*(n-1));
        }
        x = gsl_ran_exponential(rng, m);

        if(t + x < pnode->end) {
            // coalescent event within interval
            t += x;
            for(i=0; i < pnode->nsamples; ++i)
                Gene_addToBranch(pnode->sample[i], x);

            // choose a random pair to join
            i = gsl_rng_uniform_int(rng, pnode->nsamples);
            j = gsl_rng_uniform_int(rng, pnode->nsamples-1);
            if(j >= i)
                ++j;
            if(j < i) {
                k = i;
                i = j;
                j = k;
            }
            assert(i<j);

            pnode->sample[i] = Gene_join(pnode->sample[i], pnode->sample[j]);
            checkmem(pnode->sample[i], __FILE__, __LINE__);
            --pnode->nsamples;
            if(j != pnode->nsamples) {
                pnode->sample[j] = pnode->sample[pnode->nsamples];
                pnode->sample[pnode->nsamples] = NULL;
            }
        }else{
            // no coalescent event within interval
            x = pnode->end - t;
            for(i=0; i < pnode->nsamples; ++i)
                Gene_addToBranch(pnode->sample[i], x);
            t = pnode->end;
        }
    }

    // Make sure we're at the end of the interval
    if(t < pnode->end) {
        assert(pnode->nsamples < 2);
        x = pnode->end - t;
        for(i=0; i < pnode->nsamples; ++i)
            Gene_addToBranch(pnode->sample[i], x);
        t = pnode->end;
    }

    // If we have both samples and parents, then move samples to parents
    if(pnode->nsamples > 0 && pnode->nparents > 0) {
        assert(t == pnode->end);
        for(i=0; i < pnode->nsamples; ++i) {
            // move current sample to parental population 
            if(pnode->nparents > 1 && gsl_rng_uniform(rng) < pnode->mix) {
                PopNode_addSample(pnode->parent[1], pnode->sample[i]);
            }else if(pnode->nparents > 0)
                PopNode_addSample(pnode->parent[0], pnode->sample[i]);
        }
        pnode->nsamples = 0;
    }

    PopNode_sanityCheck(pnode, __FILE__, __LINE__);
    return (pnode->nsamples == 1 ? pnode->sample[0] : NULL);
}

/// Free node but not descendants
void PopNode_free(PopNode *self) {
    free(self);
}

/// survival fuction
double survival(double t, double twoN) {
    assert(t>=0.0);
    assert(twoN>0.0);
    return exp(-t/twoN);
}
