/**
 * @file popnode.c
 * @author Alan R. Rogers
 * @brief A single segment of a population tree.
 *
 * PopNode objects can be linked together into a network, which models
 * bifurcation of populations and gene flow among them. Each PopNode
 * knows its size and duration. It has pointers to parents and
 * children. If it has two parents, there is also a mixing parameter,
 * which determines what fraction of the node derives from each parent.
 *
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "popnode.h"
#include "gene.h"
#include "misc.h"
#include "network.h"
#include "nodestore.h"
#include "parstore.h"
#include "error.h"
#include "dtnorm.h"
#include <stdbool.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <gsl/gsl_randist.h>

static void  PopNode_addSample(PopNode * self, Gene * gene);
static void  PopNode_free(PopNode * self);
static void  PopNode_printShallow(PopNode * self, FILE * fp);
static void  PopNode_sanityCheck(PopNode * self, const char *file, int lineno);
static void  PopNode_sanityFromLeaf(PopNode * self, const char *file, int line);
static int   PopNode_nsamples(PopNode * self);

/// Check for errors in PopNode tree. Call this from each leaf node.
static void PopNode_sanityFromLeaf(PopNode * self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self != NULL, file, line);
    switch (self->nparents) {
    case 0:
        REQUIRE(self->parent[0] == NULL, file, line);
        REQUIRE(self->parent[1] == NULL, file, line);
        REQUIRE(self->mix == 0.0, file, line);
        REQUIRE(isinf(self->end) && self->end > 0, file, line);
        break;
    case 1:
        REQUIRE(self->parent[0] != NULL, file, line);
        REQUIRE(self->parent[1] == NULL, file, line);
        REQUIRE(self->mix == 0.0, file, line);
        REQUIRE(isfinite(self->end) && self->end >= 0, file line);
        REQUIRE(self->end == self->parent[0]->start, file, line);
        break;
    default:
        REQUIRE(self->nparents == 2, file, line);
        REQUIRE(self->parent[0] != NULL, file, line);
        REQUIRE(self->parent[1] != NULL, file, line);
        REQUIRE(isfinite(self->end), file, line);
        REQUIRE(self->mix >= 0.0, file, line);
        REQUIRE(self->end == self->parent[0]->start, file, line);
        REQUIRE(self->end == self->parent[1]->start, file, line);
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
    REQUIRE(self->start <= self->end, file, line);
    if(self->nparents > 0)
        PopNode_sanityFromLeaf(self->parent[0], file, line);
    if(self->nparents > 1)
        PopNode_sanityFromLeaf(self->parent[1], file, line);
#endif
}

/// Find root of population tree, starting from given node.
void *PopNode_root(void * vself) {
    PopNode *self = vself, *r0, *r1;
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
            fprintf(stderr, "%s:%s:%d: Population network has multiple roots\n",
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

/// Remove all references to samples from tree of populations.
/// Doesn't free the Gene objects, because they aren't owned by
/// PopNode.
void PopNode_clear(PopNode * self) {
    int         i;
    for(i = 0; i < self->nchildren; ++i)
        PopNode_clear(self->child[i]);

    self->nsamples = 0;
    //memset(self->sample, 0, sizeof(self->sample));
    PopNode_sanityCheck(self, __FILE__, __LINE__);
}

/// Return 1 if PopNode tree is empty of samples
int PopNode_isClear(const PopNode * self) {
    if(self == NULL)
        return 1;
    if(self->nsamples > 0)
        return 0;

    for(int i = 0; i < self->nchildren; ++i) {
        if(!PopNode_isClear(self->child[i]))
            return 0;
    }
    return 1;
}

/// Print a PopNode and (recursively) its descendants.
void PopNode_print(FILE * fp, void * vself, int indent) {
    PopNode *self = vself;
    for(int i = 0; i < indent; ++i)
        fputs("   ", fp);
    fprintf(fp, "%p twoN=%lf ntrval=(%lf,", self, self->twoN, self->start);
    fprintf(fp, "%lf)\n", self->end);

    for(int i = 0; i < self->nchildren; ++i)
        PopNode_print(fp, self->child[i], indent + 1);
}

/// Print a PopNode but not its descendants.
static void PopNode_printShallow(PopNode * self, FILE * fp) {
    fprintf(fp, "%p twoN=%lf ntrval=(%lf,", self, self->twoN, self->start);
    fprintf(fp, "%lf)", self->end);
    if(self->mix > 0.0)
        fprintf(fp, " mix=%lf", self->mix);

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

/// Return the number of samples in a PopNode
static int PopNode_nsamples(PopNode * self) {
    return self->nsamples;
}

/// PopNode constructor
void *PopNode_new(int twoN_i, int start_i, ParStore *ps, NodeStore * ns) {
    PopNode    *self = NodeStore_alloc(ns);
    CHECKMEM(self);

    memset(self, 0, sizeof(*self));
    self->twoN_i = twoN_i;
    self->start_i = start_i;
    self->end_i = -1;
    self->mix_i = -1;

    self->twoN = ParStore_getVal(ps, twoN_i);
    self->start = ParStore_getVal(ps, start_i);
    self->end = INFINITY;
    self->mix = 0.0;

    PopNode_sanityCheck(self, __FILE__, __LINE__);
    return self;
}

void PopNode_update(PopNode *self, ParStore *ps) {
    assert(self);
    self->twoN = ParStore_getVal(ps, self->twoN_i);
    self->start = ParStore_getVal(ps, self->start_i);
    if(self->end_i >= 0)
        self->end = ParStore_getVal(ps, self->end_i);
    if(self->mix_i >= 0)
        self->mix = ParStore_getVal(ps, self->mix_i);
    if(self->nchildren > 0)
        PopNode_updateValues(self->child[0], ps);
    if(self->nchildren > 1)
        PopNode_updateValues(self->child[1], ps);
}

/// Connect parent and child
int PopNode_addChild(void * vparent, void * vchild) {
    PopNode *parent = vparent;
    PopNode *child = vchild;
    if(parent->nchildren > 1) {
        fprintf(stderr,
                "%s:%s:%d: Can't add child because parent already has %d.\n",
                __FILE__, __func__, __LINE__, parent->nchildren);
        return TOO_MANY_CHILDREN;
    }
    if(child->nparents > 1) {
        fprintf(stderr,
                "%s:%s:%d: Can't add parent because child already has %d.\n",
                __FILE__, __func__, __LINE__, child->nparents);
        return TOO_MANY_PARENTS;
    }
    if(child->start > parent->start) {
        fprintf(stderr,
                "%s:%s:%d: Child start (%lf) must be <= parent start (%lf)\n",
                __FILE__, __func__, __LINE__, child->start, parent->start);
        return DATE_MISMATCH;
    }
    if(child->end_i == -1) {
        child->end_i = parent->start_i;
        child->end = parent->start;
    } else if(child->end_i != parent->start_i) {
            fprintf(stderr, "%s:%s:%d: Date mismatch."
                    " child->end=%p != %p = parent->start\n",
                    __FILE__, __func__, __LINE__, child->end, parent->start);
            return DATE_MISMATCH;
    }

    parent->child[parent->nchildren] = child;
    child->parent[child->nparents] = parent;
    ++parent->nchildren;
    ++child->nparents;
    PopNode_sanityCheck(parent, __FILE__, __LINE__);
    PopNode_sanityCheck(child, __FILE__, __LINE__);
    return 0;
}

/// Check sanity of PopNode
static void PopNode_sanityCheck(PopNode * self, const char *file, int lineno) {
#ifndef NDEBUG
    int         i;

    REQUIRE(self != NULL, file, lineno);

    for(i = 0; i < self->nsamples; ++i)
        REQUIRE(self->sample[i] != NULL, file, lineno);
#endif
}

/// Add a sample to a PopNode
static void PopNode_addSample(PopNode * self, Gene * gene) {
    assert(self != NULL);
    assert(gene != NULL);
    if(self->nsamples == MAXSAMP) {
        fprintf(stderr, "%s:%s:%d: Too many samples\n",
                __FILE__, __func__, __LINE__);
        exit(1);
    }
    self->sample[self->nsamples] = gene;
    ++self->nsamples;
    PopNode_sanityCheck(self, __FILE__, __LINE__);
}

/// Connect a child PopNode to two parents.
/// @param[inout] child pointer to the child PopNode
/// @param[in] mPtr pointer to the gene flow variable
/// @param[inout] introgressor pointer to the introgressing parent
/// @param[inout] native pointer to the native parent
int PopNode_mix(void * vchild, int mix_i, void * vintrogressor,
                void * vnative, ParStore *ps) {
    PopNode *child = vchild, *introgressor = vintrogressor,
        *native = vnative;

    if(introgressor->nchildren > 1) {
        fprintf(stderr,"%s:%s:%d:"
                " Can't add child because introgressor already has %d.\n",
                __FILE__, __func__, __LINE__, introgressor->nchildren);
        return TOO_MANY_CHILDREN;
    }
    if(native->nchildren > 1) {
        fprintf(stderr,"%s:%s:%d:"
                " Can't add child because native parent already has %d.\n",
                __FILE__, __func__, __LINE__, native->nchildren);
        return TOO_MANY_CHILDREN;
    }
    if(child->nparents > 0) {
        fprintf(stderr, "%s:%s:%d:"
                " Can't add 2 parents because child already has %d.\n",
                __FILE__, __func__, __LINE__, child->nparents);
        return TOO_MANY_PARENTS;
    }
    if(child->end_i >= 0) {
        if(child->end_i != introgressor->start_i) {
            fprintf(stderr,"%s:%s:%d: Date mismatch."
                    " child->end_i=%d != %d=introgressor->start_i\n",
                    __FILE__, __func__, __LINE__,
                    child->end_i, introgressor->start_i);
            return DATE_MISMATCH;
        }
        if(child->end_i != native->start_i) {
            fprintf(stderr, "%s:%s:%d: Date mismatch."
                    " child->end_i=%d != %d=native->start_i\n",
                    __FILE__, __func__, __LINE__, child->end_i,
                    native->start_i);
            return DATE_MISMATCH;
        }
    } else if(native->start_i != introgressor->start_i) {
        fprintf(stderr, "%s:%s:%d: Date mismatch."
                "native->start_i=%d != %d=introgressor->start_i\n",
                __FILE__, __func__, __LINE__,
                native->start_i, introgressor->start_i);
        return DATE_MISMATCH;
    } else {
        child->end_i = native->start_i;
        child->end = native->start;
    }

    child->parent[0] = native;
    child->parent[1] = introgressor;
    child->nparents = 2;
    child->mix_i = mix_i;
    child->mix = ParStore_getVal(ps, mix_i);
    introgressor->child[introgressor->nchildren] = child;
    ++introgressor->nchildren;
    native->child[native->nchildren] = child;
    ++native->nchildren;
    PopNode_sanityCheck(child, __FILE__, __LINE__);
    PopNode_sanityCheck(introgressor, __FILE__, __LINE__);
    PopNode_sanityCheck(native, __FILE__, __LINE__);
    return 0;
}

/// Allocates a new Gene and puts it into the array within
/// PopNode. The gene isn't owned by PopNode, however. It will
/// eventually be freed by a recursive call to Gene_free, which will
/// free the root Gene and all descendants.
void PopNode_newGene(PopNode * self, unsigned ndx) {
    assert(1 + self->nsamples < MAXSAMP);
    assert(ndx < 8 * sizeof(tipId_t));

    static const tipId_t one = 1;
    Gene       *gene = Gene_new(one << ndx);
    CHECKMEM(gene);
    self->sample[self->nsamples] = gene;
    ++self->nsamples;
    PopNode_sanityCheck(self, __FILE__, __LINE__);
}

/// Coalesce gene tree within population tree.
Gene       *PopNode_coalesce(PopNode * self, gsl_rng * rng) {

    // Because this is a network rather than a tree, we may visit the
    // same node twice. On the 2nd visit, nsamples==0, so we return
    // immediately.
    if(self->nsamples == 0)
        return NULL;

    // Coalesce children first, so that the coalescent process in
    // the current node begins with samples "inherited" from
    // children. 
    if(self->child[0])
        (void) PopNode_coalesce(self->child[0], rng);
    if(self->child[1])
        (void) PopNode_coalesce(self->child[1], rng);

    unsigned long i, j, k;
    double      x;
    double      end = self->end;
    double      t = self->start;

#ifndef NDEBUG
    if(t > end) {
        fflush(stdout);
        fprintf(stderr, "ERROR:%s:%s:%d: start=%lf > %lf=end\n",
                __FILE__, __func__, __LINE__, t, end);
        PopNode_print(stderr, self, 0);
        exit(1);
    }
#endif

    // Coalescent loop continues until only one sample is left
    // or we reach the end of the interval.
    while(self->nsamples > 1 && t < end) {
        {
            int         n = self->nsamples;
            double      mean = 2.0 * self->twoN / (n * (n - 1));
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
            CHECKMEM(self->sample[i]);
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
        x = end - t;            // may be infinite
        for(i = 0; i < self->nsamples; ++i)
            Gene_addToBranch(self->sample[i], x);
        t = end;                // may be infinite
    }
    // If we have both samples and parents, then move samples to parents
    if(self->nsamples > 0 && self->nparents > 0) {
        assert(t == end);
        assert(-1 != self->mix_i || self->nparents <= 1);
        switch (self->nparents) {
        case 1:
            // add all samples to parent 0
            for(i = 0; i < self->nsamples; ++i) {
                assert(self->sample[i]);
                PopNode_addSample(self->parent[0], self->sample[i]);
            }
            break;
        default:
            // distribute samples among parents
            assert(self->nparents == 2);
            for(i = 0; i < self->nsamples; ++i) {
                if(gsl_rng_uniform(rng) < self->mix) {
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
static void PopNode_free(PopNode * self) {
    free(self);
}

/// Return 1 if parameters satisfy inequality constraints, or 0 otherwise.
int PopNode_feasible(const PopNode * self, Bounds bnd, int verbose) {
    if(self->twoN < bnd.lo_twoN || self->twoN > bnd.hi_twoN) {
        if(verbose)
            fprintf(stderr, "%s FAIL: twoN=%lg not in [%lg, %lg]\n",
                    __func__, self->twoN, bnd.lo_twoN, bnd.hi_twoN);
        return 0;
    }

    if(self->start > bnd.hi_t || self->start < bnd.lo_t) {
        if(verbose)
            fprintf(stderr, "%s FAIL: start=%lg not in [%lg, %lg]\n",
                    __func__, self->start, bnd.lo_t, bnd.hi_t);
        return 0;
    }

    switch (self->nparents) {
    case 2:
        if(self->start > self->parent[1]->start) {
            if(verbose)
                fprintf(stderr, "%s FAIL: child=%lg older than parent=%lg\n",
                        __func__, self->start, self->parent[1]->start);
            return 0;
        }
        // fall through
    case 1:
        if(self->start > self->parent[0]->start) {
            if(verbose)
                fprintf(stderr, "%s FAIL: child=%lg older than parent=%lg\n",
                        __func__, self->start, self->parent[0]->start);
            return 0;
        }
        break;
    default:
        break;
    }

    switch (self->nchildren) {
    case 2:
        if(self->start < self->child[1]->start) {
            if(verbose)
                fprintf(stderr,
                        "%s FAIL: parent=%lg younger than child=%lg\n",
                        __func__, self->start, self->child[1]->start);
            return 0;
        }
        // fall through
    case 1:
        if(self->start < self->child[0]->start) {
            if(verbose)
                fprintf(stderr,
                        "%s FAIL: parent=%lg younger than child=%lg\n",
                        __func__, self->start, self->child[0]->start);
            return 0;
        }
        break;
    default:
        break;
    }

    if(self->mix_i != -1) {
        if(self->mix < 0.0 || self->mix > 1.0) {
            if(verbose)
                fprintf(stderr, "%s FAIL: mix=%lg not in [0, 1]\n",
                        __func__, self->mix);
            return 0;
        }
    }

    for(int i = 0; i < self->nchildren; ++i) {
        if(0 == PopNode_feasible(self->child[i], bnd, verbose))
            return 0;
    }
    return 1;
}

/// Add dp to each parameter pointer, using ordinary (not pointer)
/// arithmetic.
void PopNode_shiftParamPtrs(PopNode * self, size_t dp, int sign) {
    SHIFT_PTR(self->twoN, dp, sign);
    SHIFT_PTR(self->start, dp, sign);
    SHIFT_PTR(self->end, dp, sign);
    SHIFT_PTR(self->mix, dp, sign);
}

/// Add dp to each PopNode pointer, using ordinary (not pointer)
/// arithmetic.
void PopNode_shiftPopNodePtrs(PopNode * self, size_t dp, int sign) {
    int         i;
    for(i = 0; i < self->nparents; ++i)
        SHIFT_PTR(self->parent[i], dp, sign);

    for(i = 0; i < self->nchildren; ++i)
        SHIFT_PTR(self->child[i], dp, sign);
}

#ifdef TEST

#  include <string.h>
#  include <assert.h>
#  include <time.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

int main(int argc, char **argv) {

    int         verbose = 0;

    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xpopnode [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    Network_init(SIM);

    int         nseg = 10;
    PopNode     v[nseg];
    NodeStore  *ns = NodeStore_new(nseg, sizeof(v[0]), v);
    CHECKMEM(ns);

    double      twoN0 = 1.0, start0 = 0.0;
    PopNode    *p0 = PopNode_new(&twoN0, &start0, ns);
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

    double      twoN1 = 100.0, start1 = 123.0;
    PopNode    *p1 = PopNode_new(&twoN1, &start1, ns);
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

    tipId_t     id1 = 1;
    tipId_t     id2 = 2;
    Gene       *g1 = Gene_new(id1);
    Gene       *g2 = Gene_new(id2);
    PopNode_addSample(p1, g1);
    PopNode_addSample(p1, g2);
    assert(p1->nsamples == 2);

    int status=PopNode_addChild(p1, p0);
    assert(status==0);
    assert(p1->nchildren == 1);
    assert(p0->nparents == 1);
    assert(p1->child[0] == p0);
    assert(p0->parent[0] == p1);

    if(verbose) {
        PopNode_printShallow(p1, stdout);
        PopNode_printShallow(p0, stdout);
    }

    ParStore   *ps = ParStore_new();

    size_t      twoNloc = (size_t) p1->twoN;
    size_t      startloc = (size_t) p1->start;
    size_t      endloc = (size_t) p1->end;
    PopNode_shiftParamPtrs(p1, (size_t) 1u, 1);
    assert(endloc == 0u || twoNloc + 1u == (size_t) p1->twoN);
    assert(endloc == 0u || startloc + 1u == (size_t) p1->start);
    assert(endloc == 0u || endloc + 1u == (size_t) p1->end);

    int         i;
    size_t      parent[2], child[2];
    for(i = 0; i < p1->nparents; ++i)
        parent[i] = (size_t) p1->parent[i];
    for(i = 0; i < p1->nchildren; ++i)
        child[i] = (size_t) p1->child[i];
    PopNode_shiftPopNodePtrs(p1, (size_t) 1u, 1);
    for(i = 0; i < p1->nparents; ++i)
        assert(parent[i] + 1u == (size_t) p1->parent[i]);
    for(i = 0; i < p1->nchildren; ++i)
        assert(child[i] + 1u == (size_t) p1->child[i]);

    unitTstResult("PopNode", "untested");

    NodeStore_free(ns);

    ParStore_free(ps);
    Gene_free(g1);
    Gene_free(g2);

    return 0;
}
#endif
