/**
 * @file gptree.c
 * @brief Methods for simulating gene genealogies within a given tree
 * of populations, and allowing populations to mix and also to split.
 */

#include "gptree.h"
#include "misc.h"
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
    double K;          /* ratio of current pop size to ancestral pop size */
    double start, end; /* duration of this PopNode */
    struct PopNode *parent[2];
    struct PopNode *child[2];
    double mix;  /*mix=frac of pop derived from parent[1] */
    Gene *sample[MAXSAMPLES];
};

/// Find root of population tree, starting from given node.
PopNode *PopNode_root(PopNode *self) {
    assert(self);
    switch(self->nparents) {
    case 0:
        return self;
        break;
    case 1:
        return PopNode_root(self->parent[0]);
        break;
    case 2:
        PopNode *r0 = PopNode_root(self->parent[0]);
        PopNode *r1 = PopNode_root(self->parent[1]);
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

static void PopNode_sanityCheck(PopNode *pnode, const char *file, int lineno);

void PopNode_set(PopNode *self, double K, double start, double end) {
    self->K = K;
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
        fputs("---|", fp);
    fprintf(fp, "%p K=%lf ntrval=(%lf,",
            pnode, pnode->K, pnode->start);
    if(pnode->end < DBL_MAX)
        fprintf(fp, "%lf)\n", pnode->end);
    else
        fprintf(fp, "Inf)\n");

    for(i=0; i < pnode->nchildren; ++i)
        PopNode_print(fp, pnode->child[i], indent+1);
}

Gene *Gene_new(tipId_t tipId) {
    Gene *gene = malloc(sizeof(Gene));
    checkmem(gene, __FILE__, __LINE__);

    gene->tipId = tipId;
    gene->parent = gene->lchild = gene->rchild = NULL;
    gene->branch = 0.0;

    return gene;
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
    return parent;
}

/**
 * Return the length of the deepest coalescent interval.
 * Set *tipId equal to the tipId of the oldest of the two children
 * of the given node.
 */
double Gene_lastInterval(Gene *gene, tipId_t *tipId) {
    if(gene->lchild==NULL || gene->rchild==NULL)
        eprintf("%s:%s%d: root has null children: can't calculate branch length\n",
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

PopNode *PopNode_new(double K, double start, double end) {
    PopNode *pnode = malloc(sizeof(PopNode));
    checkmem(pnode, __FILE__, __LINE__);

    pnode->nparents = pnode->nchildren = pnode->nsamples = 0;
    pnode->K = K;
    pnode->start = start;
    pnode->end = end;

    memset(pnode->sample, 0, sizeof(pnode->sample));
    memset(pnode->parent, 0, sizeof(pnode->parent));
    memset(pnode->child, 0, sizeof(pnode->child));

    PopNode_sanityCheck(pnode, __FILE__, __LINE__);
    return pnode;
}

/** Add a child population to parent population */
void PopNode_addChild(PopNode *parent, PopNode *child) {
    assert(parent->nchildren < 2);
    parent->child[parent->nchildren] = child;
    ++parent->nchildren;
    ++child->nparents;
    assert(child->nparents < 3);
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
    if(pnode->nsamples == MAXSAMPLES) {
        fprintf(stderr,"%s:%s:%d: Too many samples\n",
                __FILE__,__func__,__LINE__);
        exit(1);
    }
    pnode->sample[pnode->nsamples] = gene;
    ++pnode->nsamples;
    PopNode_sanityCheck(pnode, __FILE__, __LINE__);
}

void PopNode_mix(PopNode *pnode, double m, PopNode *immigrant, PopNode *native) {
    pnode->parent[0] = native;
    pnode->parent[1] = immigrant;
    pnode->mix = m;
    PopNode_addChild(immigrant, pnode);
    PopNode_addChild(native, pnode);
    assert(pnode->nparents == 2);
    PopNode_sanityCheck(pnode, __FILE__, __LINE__);
    PopNode_sanityCheck(immigrant, __FILE__, __LINE__);
    PopNode_sanityCheck(native, __FILE__, __LINE__);
}

void PopNode_endToEnd(PopNode *pnode, PopNode *ancestor) {
    pnode->nparents = 1;
    pnode->parent[0] = ancestor;
    pnode->mix = 0.0;
    PopNode_addChild(ancestor, pnode);
    PopNode_sanityCheck(pnode, __FILE__, __LINE__);
    PopNode_sanityCheck(ancestor, __FILE__, __LINE__);
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
    assert(1 + pnode->nsamples < MAXSAMPLES);

    assert(ndx < sizeof(tipId_t));
    Gene *gene = Gene_new(1UL << ndx);
    checkmem(gene, __FILE__, __LINE__);
    pnode->sample[pnode->nsamples] = gene;
    ++pnode->nsamples;
    PopNode_sanityCheck(pnode, __FILE__, __LINE__);
}

/**
 * Coalesce gene tree within population tree.
 */
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

    /*
     * Coalescent loop continues until only one sample is left
     * or we reach the end of the interval.
     */
    while(pnode->nsamples > 1 && t < pnode->end) {
        double m;
        {
            int n = pnode->nsamples;
            m = 2.0 * pnode->K / (n*(n-1));
        }
        x = gsl_ran_exponential(rng, m);

        if(t + x < pnode->end) {
            /* coalescent event within interval */
            t += x;
            for(i=0; i < pnode->nsamples; ++i)
                Gene_addToBranch(pnode->sample[i], x);
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
            /* no coalescent event within interval */
            x = pnode->end - t;
            for(i=0; i < pnode->nsamples; ++i)
                Gene_addToBranch(pnode->sample[i], x);
            t = pnode->end;
        }
    }

    /* Make sure we're at the end of the interval */
    if(t < pnode->end) {
        assert(pnode->nsamples < 2);
        x = pnode->end - t;
        for(i=0; i < pnode->nsamples; ++i)
            Gene_addToBranch(pnode->sample[i], x);
        t = pnode->end;
    }

    /* If we have both samples and parents, then move samples to parents */
    if(pnode->nsamples > 0 && pnode->nparents > 0) {
        assert(t == pnode->end);
        for(i=0; i < pnode->nsamples; ++i) {
            /* move current sample to parental population */
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

/** Need to figure out how to free gene tree and pop tree. */
void PopNode_free(PopNode *pnode) {
    PopNode_sanityCheck(pnode, __FILE__, __LINE__);
    free(pnode);
}

/** survival fuction */
double survival(double t, double K) {
    assert(t>=0.0);
    assert(K>0.0);
    return exp(-t/K);
}
