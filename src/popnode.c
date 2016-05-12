/**
 * @file popnode.c
 * @brief Class PopNode represents a single segment of a population
 * tree. These nodes are linked together into a network, which models
 * bifurcation of populations and gene flow among them.
 */

#include "popnode.h"
#include "gene.h"
#include "misc.h"
#include <string.h>
#include <gsl/gsl_randist.h>

struct NodeStore {
    int nused, len;
    PopNode *v; // not locally owned
};

static void PopNode_sanityCheck(PopNode * self, const char *file, int lineno);

void PopNode_sanityFromLeaf(PopNode * self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self != NULL, file, line);
    switch (self->nparents) {
    case 0:
        REQUIRE(self->parent[0] == NULL, file, line);
        REQUIRE(self->parent[1] == NULL, file, line);
        REQUIRE(self->mix == NULL, file, line);
        REQUIRE(self->end == NULL, file, line);
        break;
    case 1:
        REQUIRE(self->parent[0] != NULL, file, line);
        REQUIRE(self->parent[1] == NULL, file, line);
        REQUIRE(self->mix == NULL, file, line);
        break;
    default:
        REQUIRE(self->nparents == 2, file, line);
        REQUIRE(self->parent[0] != NULL, file, line);
        REQUIRE(self->parent[1] != NULL, file, line);
		REQUIRE(self->end != NULL, file, line);
		REQUIRE(self->mix != NULL, file, line);
        REQUIRE(*self->mix >= 0.0, file, line);
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
    REQUIRE(self->end==NULL || *self->start <= *self->end, file, line);
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

/// Return 1 if PopNode tree is empty of samples
int PopNode_isClear(const PopNode *self) {
    if(self==0)
        return 1;
    if(self->nsamples > 0)
        return 0;

    int i;
    for(i=0; i < self->nchildren; ++i) {
        if(!PopNode_isClear(self->child[i]))
            return 0;
    }
    return 1;
}

void PopNode_print(FILE * fp, PopNode * self, int indent) {
    int         i;
    for(i = 0; i < indent; ++i)
        fputs("   ", fp);
    fprintf(fp, "%p twoN=%lf ntrval=(%lf,", self, *self->twoN,
			*self->start);
    if(self->end != NULL)
        fprintf(fp, "%lf)\n", *self->end);
    else
        fprintf(fp, "Inf)\n");

    for(i = 0; i < self->nchildren; ++i)
        PopNode_print(fp, self->child[i], indent + 1);
}

void PopNode_printShallow(PopNode * self, FILE * fp) {
    fprintf(fp, "%p twoN=%lf ntrval=(%lf,",
            self, *self->twoN, *self->start);
    if(self->end != NULL)
        fprintf(fp, "%lf)", *self->end);
    else
        fprintf(fp, "Inf)");
	if(self->mix != NULL)
		fprintf(fp, " mix=%lf", *self->mix);

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

PopNode    *PopNode_new(double *twoN, double *start, NodeStore *ns) {
    PopNode    *new = NodeStore_alloc(ns);
    checkmem(new, __FILE__, __LINE__);

    new->nparents = new->nchildren = new->nsamples = 0;
    new->twoN = twoN;
    new->mix = NULL;
    new->start = start;
    new->end = NULL;

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
    if(*child->start > *parent->start)
        eprintf("%s:%s:%d: Child start (%lf) must be <= parent start (%lf)\n",
                __FILE__, __func__, __LINE__, *child->start, *parent->start);
    if(child->end == NULL) {
        child->end = parent->start;
    } else {
		if(child->end != parent->start)
			eprintf("%s:%s:%d: Date mismatch."
					" child->end=%p != %p = parent->start\n",
					__FILE__, __func__, __LINE__,
					child->end, parent->start);
	}
    parent->child[parent->nchildren] = child;
    child->parent[child->nparents] = parent;
    ++parent->nchildren;
    ++child->nparents;
    PopNode_sanityCheck(parent, __FILE__, __LINE__);
    PopNode_sanityCheck(child, __FILE__, __LINE__);
}

static void PopNode_sanityCheck(PopNode * self, const char *file,
								int lineno) {
#ifndef NDEBUG
    int         i;

    REQUIRE(self != NULL, file, lineno);

    for(i = 0; i < self->nsamples; ++i)
        REQUIRE(self->sample[i] != NULL, file, lineno);
#endif
}

void PopNode_addSample(PopNode * self, Gene * gene) {
	assert(self!=NULL);
	assert(gene!=NULL);
    if(self->nsamples == MAXSAMP) {
        fprintf(stderr, "%s:%s:%d: Too many samples\n",
                __FILE__, __func__, __LINE__);
        exit(1);
    }
    self->sample[self->nsamples] = gene;
    ++self->nsamples;
    PopNode_sanityCheck(self, __FILE__, __LINE__);
}

void PopNode_mix(PopNode * child, double *mPtr, PopNode * introgressor,
                 PopNode * native) {
    if(introgressor->nchildren > 1)
        eprintf("%s:%s:%d:"
				" Can't add child because introgressor already has %d.\n",
				 __FILE__, __func__, __LINE__, introgressor->nchildren);  
    if(native->nchildren > 1)
        eprintf("%s:%s:%d:"
				" Can't add child because native parent already has %d.\n",
				 __FILE__, __func__, __LINE__, native->nchildren);
    if(child->nparents > 0)
        eprintf("%s:%s:%d:"
				" Can't add 2 parents because child already has %d.\n",
				 __FILE__, __func__, __LINE__, child->nparents);
    if(child->end != NULL) {
        if(child->end != introgressor->start)
            eprintf("%s:%s:%d: Date mismatch."
					 " child->end=%p != %p=introgressor->start\n",
					 __FILE__, __func__, __LINE__,
					child->end, introgressor->start);
        if(child->end != native->start)
            eprintf("%s:%s:%d: Date mismatch."
					 " child->end=%p != %p=native->start\n",
					 __FILE__, __func__, __LINE__,
					child->end, native->start);
	} else if(native->start != introgressor->start) {
		eprintf("%s:%s:%d: Date mismatch."
				 "native->start=%p != %p=introgressor->start\n",
				 __FILE__, __func__, __LINE__,
				 native->start, introgressor->start);
    } else
        child->end = native->start;

    child->parent[0] = native;
    child->parent[1] = introgressor;
    child->nparents = 2;
    child->mix = mPtr;
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

    assert(ndx < 8*sizeof(tipId_t));
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
	double end = (NULL==self->end ? HUGE_VAL : *self->end);

    if(self->child[0])
        (void) PopNode_coalesce(self->child[0], rng);
    if(self->child[1])
        (void) PopNode_coalesce(self->child[1], rng);

    double      t = *self->start;
#ifndef NDEBUG
    if(t > end) {
        PopNode_print(stdout, self, 0);
		eprintf("%s:%s:%d: start=%lf > %lf=end\n",
				__FILE__,__func__,__LINE__, t, end);
	}
#endif

    // Coalescent loop continues until only one sample is left
    // or we reach the end of the interval.
    while(self->nsamples > 1 && t < end) {
        {
            int         n = self->nsamples;
            double      mean = 2.0 * *self->twoN / (n * (n - 1));
			x = gsl_ran_exponential(rng, mean);
        }

        if(t + x < end) {
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
			assert(isfinite(end));
            x = end - t;
            for(i = 0; i < self->nsamples; ++i)
                Gene_addToBranch(self->sample[i], x);
            t = end;
        }
    }

    // Make sure we're at the end of the interval
    if(t < end) {
        assert(self->nsamples < 2);
        x = end - t;  // may be infinite
        for(i = 0; i < self->nsamples; ++i)
            Gene_addToBranch(self->sample[i], x);
        t = end;      // may be infinite
    }

    // If we have both samples and parents, then move samples to parents
    if(self->nsamples > 0 && self->nparents > 0) {
        assert(t == end);
		assert(NULL!=self->mix || self->nparents <= 1);
		switch(self->nparents) {
		case 1:
			// add all samples to parent 0
			for(i = 0; i < self->nsamples; ++i) {
				assert(self->sample[i]);
                PopNode_addSample(self->parent[0], self->sample[i]);
			}
			break;
		default:
			// distribute samples among parents
			assert(self->nparents==2);
			for(i = 0; i < self->nsamples; ++i) {
				if(gsl_rng_uniform(rng) < *self->mix) {
					assert(self->sample[i]);
					PopNode_addSample(self->parent[1], self->sample[i]);
				} else {
					assert(self->sample[i]);
					PopNode_addSample(self->parent[0], self->sample[i]);
				}
			}
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

/// Add dp to each parameter pointer, using ordinary (not pointer)
/// arithmetic.
void PopNode_shiftParamPtrs(PopNode *self, size_t dp) {
    INCR_PTR(self->twoN, dp);
    INCR_PTR(self->start, dp);
    INCR_PTR(self->end, dp);
}

/// Add dp to each PopNode pointer, using ordinary (not pointer)
/// arithmetic.
void PopNode_shiftPopNodePtrs(PopNode *self, size_t dp) {
    int i;
    for(i=0; i < self->nparents; ++i)
        INCR_PTR(self->parent[i], dp);

    for(i=0; i < self->nchildren; ++i)
        INCR_PTR(self->child[i], dp);
}

NodeStore *NodeStore_new(int len, PopNode *v) {
    NodeStore *self = malloc(sizeof(NodeStore));
    CHECKMEM(self);

    self->nused = 0;
    self->len = len;
    self->v = v;
    return self;
}

void NodeStore_free(NodeStore *self) {
    // Does not free self->v
    free(self);
}

PopNode *NodeStore_alloc(NodeStore *self) {
    if(self->nused >= self->len)
        eprintf("%s:%s:%d: Ran out of PopNode objects.\n",
                __FILE__, __func__, __LINE__);
    return &self->v[self->nused++];
}

/// Set everything to zero.
void SampNdx_init(SampNdx * self) {
    memset(self, 0, sizeof(*self));
}

/// Add samples for a single population. Should be called once for
/// each sampled population.
void SampNdx_addSamples(SampNdx * self, unsigned nsamples,
						PopNode * pnode) {
    unsigned    i;
    if(self->n + nsamples >= MAXSAMP)
        eprintf("%s:%s:%d: too many samples\n", __FILE__, __func__, __LINE__);
    for(i = 0; i < nsamples; ++i) {
        self->node[self->n] = pnode;
        self->n += 1;
    }
}

/// Put samples into the gene tree. Should be done at the start of
/// each simulation.
void SampNdx_populateTree(SampNdx * self) {
    unsigned    i;
    for(i = 0; i < self->n; ++i)
        PopNode_newGene(self->node[i], i);
}

unsigned SampNdx_size(SampNdx * self) {
    return self->n;
}

/// This equality check doesn't do much, because the pointers
/// in different SampNdx objects don't have to be (in fact shouldn't
/// be) equal.
int         SampNdx_equals(SampNdx *lhs, SampNdx *rhs){
    if(lhs==NULL && rhs==NULL)
        return 1;
    if(lhs==NULL || rhs==NULL)
        return 0;
    if(lhs->n != rhs->n)
        return 0;
    return 1;
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
            fprintf(stderr,"usage: xpopnode [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    tipId_t id1 = 1;
    tipId_t id2 = 2;

    int nseg = 10;
    PopNode v[nseg];
    NodeStore *ns = NodeStore_new(nseg, v);
    CHECKMEM(ns);

    PopNode *node = NodeStore_alloc(ns);
    assert(node == v);
    node = NodeStore_alloc(ns);
    assert(node == &v[1]);
    node = NodeStore_alloc(ns);
    assert(node == &v[2]);

    NodeStore_free(ns);
	unitTstResult("NodeStore", "OK");

    ns = NodeStore_new(nseg, v);
    CHECKMEM(ns);

    double twoN0 = 1.0, start0= 0.0;
    PopNode *p0 = PopNode_new(&twoN0, &start0, ns);
    assert(p0->twoN == &twoN0);
    assert(p0->start == &start0);
    assert(p0->end == NULL);
    assert(p0->mix == NULL);
    assert(p0->nsamples == 0);
    assert(p0->nchildren == 0);
    assert(p0->child[0] == NULL);
    assert(p0->child[1] == NULL);
    assert(p0->parent[0] == NULL);
    assert(p0->parent[1] == NULL);

	if(verbose) 
		PopNode_printShallow(p0, stdout);

    double twoN1 = 100.0, start1= 123.0;
    PopNode *p1 = PopNode_new(&twoN1, &start1, ns);
    assert(p1->twoN == &twoN1);
    assert(p1->start == &start1);
    assert(p1->end == NULL);
    assert(p1->mix == NULL);
    assert(p1->nsamples == 0);
    assert(p1->nchildren == 0);
    assert(p1->child[0] == NULL);
    assert(p1->child[1] == NULL);
    assert(p1->parent[0] == NULL);
    assert(p1->parent[1] == NULL);

    Gene *g1 = Gene_new(id1);
    Gene *g2 = Gene_new(id2);
    PopNode_addSample(p1, g1);
    PopNode_addSample(p1, g2);
    assert(p1->nsamples == 2);

    PopNode_addChild(p1, p0);
    assert(p1->nchildren == 1);
    assert(p0->nparents == 1);
    assert(p1->child[0] == p0);
    assert(p0->parent[0] == p1);

	if(verbose) 
		PopNode_printShallow(p1, stdout);

    size_t twoNloc = (size_t) p1->twoN;
    size_t startloc = (size_t) p1->start;
    size_t endloc = (size_t) p1->end;
    PopNode_shiftParamPtrs(p1, (size_t) 1u);
    assert(twoNloc+1u == (size_t) p1->twoN);
    assert(startloc+1u == (size_t) p1->start);
    assert(endloc+1u == (size_t) p1->end);

    int i;
    size_t parent[2], child[2];
    for(i=0; i < p1->nparents; ++i)
        parent[i] = (size_t) p1->parent[i];
    for(i=0; i < p1->nchildren; ++i)
        child[i] = (size_t) p1->child[i];
    PopNode_shiftPopNodePtrs(p1, (size_t) 1u);
    for(i=0; i < p1->nparents; ++i)
        assert(parent[i]+1u == (size_t) p1->parent[i]);
    for(i=0; i < p1->nchildren; ++i)
        assert(child[i]+1u == (size_t) p1->child[i]);

    unitTstResult("PopNode", "untested");

    SampNdx     sndx = {.n = 3 };
    assert(sndx.n == 3);
    assert(SampNdx_size(&sndx) == 3);

    SampNdx_init(&sndx);
    assert(SampNdx_size(&sndx) == 0);

	double twoN = 100.0;
	double start = 20.0;
    PopNode    *pnode = PopNode_new(&twoN, &start, ns);
    SampNdx_addSamples(&sndx, 1, pnode);
    SampNdx_addSamples(&sndx, 2, pnode);

    assert(3 == SampNdx_size(&sndx));
    SampNdx_populateTree(&sndx);
    assert(3 == PopNode_nsamples(pnode));
    NodeStore_free(ns);

    SampNdx     sndx2 = {.n = 3 };
    SampNdx_init(&sndx2);
	double twoN2 = 100.0;
	double start2 = 20.0;
    PopNode v2[nseg];
    NodeStore *ns2 = NodeStore_new(nseg, v2);
    CHECKMEM(ns);
    pnode = PopNode_new(&twoN2, &start2, ns2);
    SampNdx_addSamples(&sndx2, 1, pnode);
    SampNdx_addSamples(&sndx2, 2, pnode);
    SampNdx_populateTree(&sndx2);
    NodeStore_free(ns2);
    assert(SampNdx_equals(&sndx, &sndx2));

	unitTstResult("SampNdx", "OK");

    return 0;
}
#endif
