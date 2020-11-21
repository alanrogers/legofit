
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
#include "parstore.h"
#include "ptrptrmap.h"
#include "error.h"
#include <stdbool.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <gsl/gsl_randist.h>

static void PopNode_transferSample(PopNode * self, Gene * gene);
static void PopNode_printShallow(PopNode * self, FILE * fp);
static void PopNode_sanityCheck(PopNode * self, const char *file, int lineno);
static void PopNode_sanityFromLeaf(PopNode * self, const char *file,
                                   int line);
static int PopNode_nsamples(PopNode * self);
static void PopNode_duplicate_nodes(PopNode * old, PtrPtrMap * ppm);
static void unlink_child(PopNode * child, PopNode * parent);
static int PopNode_equals_r(PopNode * a, PopNode * b);
static void PopNode_print_r(PopNode *self, FILE * fp, int indent);

struct PopNode {
    char *label;                // name of current segment
    int nparents, nchildren, nsamples;
    double twoN;                // haploid pop size
    double start, end;          // duration of this PopNode
    double mix;                 // frac of pop derived from parent[1]
    int visited;                // has the coalescent visited this node yet?

    // indices into ParStore array
    int twoN_i, start_i, end_i, mix_i;

    struct PopNode *parent[2];
    struct PopNode *child[2];

    Gene *sample[MAXSAMP];      // not locally owned
};

int PopNode_equals(PopNode * a, PopNode * b) {
    PopNode_unvisit(a);
    PopNode_unvisit(b);
    return PopNode_equals_r(a, b);
}

static int PopNode_equals_r(PopNode * a, PopNode * b) {
    if(a->visited != b->visited)
        return 0;
    if(a->visited) {
        // We've been here before. Had the two nodes been unequal last
        // time, the algorithm would not have visited this spot a 2nd
        // time. So because we're here, the two nodes must be equal.
        return 1;
    }
    a->visited = b->visited = 1;

    if(a->nparents != b->nparents)
        return 0;
    if(a->nchildren != b->nchildren)
        return 0;
    if(a->twoN != b->twoN)
        return 0;
    if(a->start != b->start)
        return 0;
    if(a->end != b->end)
        return 0;
    if(a->mix != b->mix)
        return 0;
    if(a->twoN_i != b->twoN_i)
        return 0;
    if(a->start_i != b->start_i)
        return 0;
    if(a->end_i != b->end_i)
        return 0;
    if(a->mix_i != b->mix_i)
        return 0;
    if(a->end != b->end)
        return 0;
    for(int i = 0; i < a->nchildren; ++i) {
        if(!PopNode_equals_r(a->child[i], b->child[i]))
            return 0;
    }
    return 1;
}

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
        REQUIRE(isfinite(self->end) && self->end >= 0, file, line);
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

/// Find root of population network, starting from given node.
void *PopNode_root(void *vself) {
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
            fprintf(stderr,
                    "%s:%s:%d: Population network has multiple roots\n",
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
/// PopNode. Sets "visited" to 0 in every node.
void PopNode_clear(PopNode * self) {
    assert(self);
    for(int i = 0; i < self->nchildren; ++i) {
        assert(self->child[i]);
        PopNode_clear(self->child[i]);
    }

    self->nsamples = 0;
    self->visited = 0;
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

/// Print a network of PopNodes.
void PopNode_print(void *vroot, FILE * fp, int indent) {
    PopNode *root = vroot;
    PopNode_unvisit(root);
    PopNode_print_r(root, fp, indent);
}

/// Print a PopNode and (recursively) its descendants.
static void PopNode_print_r(PopNode *self, FILE * fp, int indent) {

    if(self->visited)
        return;
    self->visited = 1;

    stutter('|', indent, fp);
    fprintf(fp, "%s: twoN=%lg ntrval=(%lf,", self->label, self->twoN,
            self->start);
    fprintf(fp, "%lf)\n", self->end);

    if(self->nparents) {
        stutter('|', indent, fp);
        fprintf(fp, "  parents:");
        for(int i = 0; i < self->nparents; ++i)
            fprintf(fp, " %s", self->parent[i]->label);
        putc('\n', fp);
    }

    if(self->nchildren) {
        stutter('|', indent, fp);
        fprintf(fp, "  children:");
        for(int i = 0; i < self->nchildren; ++i)
            fprintf(fp, " %s", self->child[i]->label);
        putc('\n', fp);
    }

    if(self->nsamples) {
        stutter('|', indent, fp);
        fprintf(fp, "  nsamples=%d\n", self->nsamples);
    }

    for(int i = 0; i < self->nchildren; ++i)
        PopNode_print(self->child[i], fp, indent + 1);
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

/// Set all "visited" flags to false.
void PopNode_unvisit(PopNode * self) {
    if(self->nchildren > 0)
        PopNode_unvisit(self->child[0]);
    if(self->nchildren > 1)
        PopNode_unvisit(self->child[1]);
    self->visited = 0;
}

/// PopNode constructor
void *PopNode_new(int twoN_i, int start_i, ParStore * ps, const char *label) {
    PopNode *self = malloc(sizeof(PopNode));
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
    self->label = strdup(label);
    CHECKMEM(self->label);

    PopNode_sanityCheck(self, __FILE__, __LINE__);
    return self;
}

void PopNode_update(PopNode * self, ParStore * ps) {
    assert(self);
    self->twoN = ParStore_getVal(ps, self->twoN_i);
    self->start = ParStore_getVal(ps, self->start_i);
    if(self->end_i >= 0)
        self->end = ParStore_getVal(ps, self->end_i);
    if(self->mix_i >= 0)
        self->mix = ParStore_getVal(ps, self->mix_i);
    if(self->nchildren > 0)
        PopNode_update(self->child[0], ps);
    if(self->nchildren > 1)
        PopNode_update(self->child[1], ps);
}

/// Connect parent and child
int PopNode_addChild(void *vparent, void *vchild) {
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
        fprintf(stderr, "%s:%s:%d: Date mismatch.\n"
                "  child->end_i=%d != %d = parent->start_i\n",
                __FILE__, __func__, __LINE__, child->end_i, parent->start_i);
        fprintf(stderr, "  child->end=%lg != %lg = parent->start\n",
                child->end, parent->start);
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
    int i;

    REQUIRE(self != NULL, file, lineno);

    for(i = 0; i < self->nsamples; ++i)
        REQUIRE(self->sample[i] != NULL, file, lineno);
#endif
}

/// Allocates a new Gene and puts it into the array within
/// PopNode. The gene isn't owned by PopNode, however. It will
/// eventually be freed by a recursive call to Gene_free, which will
/// free the root Gene and all descendants.
void PopNode_newSample(PopNode * self, unsigned ndx) {
    assert(1 + self->nsamples < MAXSAMP);
    assert(ndx < 8 * sizeof(tipId_t));

    static const tipId_t one = 1;
    Gene *gene = Gene_new(one << ndx);
    CHECKMEM(gene);
    self->sample[self->nsamples] = gene;
    ++self->nsamples;
    PopNode_sanityCheck(self, __FILE__, __LINE__);
}

/// Transfer an existing sample to a PopNode
static void PopNode_transferSample(PopNode * self, Gene * gene) {
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
int PopNode_mix(void *vchild, int mix_i, void *vintrogressor,
                void *vnative, ParStore * ps) {
    PopNode *child = vchild, *introgressor = vintrogressor, *native = vnative;

    if(introgressor->nchildren > 1) {
        fprintf(stderr, "%s:%s:%d:"
                " Can't add child because introgressor already has %d.\n",
                __FILE__, __func__, __LINE__, introgressor->nchildren);
        return TOO_MANY_CHILDREN;
    }
    if(native->nchildren > 1) {
        fprintf(stderr, "%s:%s:%d:"
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
            fprintf(stderr, "%s:%s:%d: Date mismatch\n"
                    "  child->end_i=%d != %d=introgressor->start_i\n",
                    __FILE__, __func__, __LINE__,
                    child->end_i, introgressor->start_i);
            fprintf(stderr, "  child->end=%lg != %lg=introgressor->start\n",
                    child->end, introgressor->start);
            return DATE_MISMATCH;
        }
        if(child->end_i != native->start_i) {
            fprintf(stderr, "%s:%s:%d: Date mismatch\n"
                    "  child->end_i=%d != %d=native->start_i\n",
                    __FILE__, __func__, __LINE__, child->end_i,
                    native->start_i);
            fprintf(stderr, "  child->end=%lg != %lg=native->start\n",
                    child->end, native->start);
            return DATE_MISMATCH;
        }
    } else if(native->start_i != introgressor->start_i) {
        fprintf(stderr, "%s:%s:%d: Date mismatch\n"
                "  native->start_i=%d != %d=introgressor->start_i\n",
                __FILE__, __func__, __LINE__,
                native->start_i, introgressor->start_i);
        fprintf(stderr, "  native->start=%lg != %lg=introgressor->start\n",
                native->start, introgressor->start);
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

/// Coalesce gene tree within population tree.
Gene *PopNode_coalesce(PopNode * self, gsl_rng * rng,
                       long unsigned *event_counter) {

    UNUSED(event_counter);

    // Early return if this node has been visited already.
    if(self->visited)
        return NULL;

    // Coalesce children first, so that the coalescent process in
    // the current node begins with samples "inherited" from
    // children. 
    if(self->child[0])
        (void) PopNode_coalesce(self->child[0], rng, event_counter);
    if(self->child[1])
        (void) PopNode_coalesce(self->child[1], rng, event_counter);

    unsigned long i, j, k;
    double x;
    double end = self->end;
    double t = self->start;

#ifndef NDEBUG
    if(t > end) {
        fflush(stdout);
        fprintf(stderr, "ERROR:%s:%s:%d: start=%lf > %lf=end\n",
                __FILE__, __func__, __LINE__, t, end);
        PopNode_unvisit(self);
        PopNode_print(self, stderr, 0);
        exit(1);
    }
#endif

    // Coalescent loop continues until only one sample is left
    // or we reach the end of the interval.
    while(self->nsamples > 1 && t < end) {
        {
            int n = self->nsamples;
            double mean = 2.0 * self->twoN / (n * (n - 1));
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
                PopNode_transferSample(self->parent[0], self->sample[i]);
            }
            break;
        default:
            // distribute samples among parents
            assert(self->nparents == 2);
            for(i = 0; i < self->nsamples; ++i) {
                if(gsl_rng_uniform(rng) < self->mix) {
                    assert(self->sample[i]);
                    PopNode_transferSample(self->parent[1], self->sample[i]);
                } else {
                    assert(self->sample[i]);
                    PopNode_transferSample(self->parent[0], self->sample[i]);
                }
            }
        }
        self->nsamples = 0;
    }

    PopNode_sanityCheck(self, __FILE__, __LINE__);
    self->visited = 1;
    return (self->nsamples == 1 ? self->sample[0] : NULL);
}

/// Remove child from parent
static void unlink_child(PopNode * child, PopNode * parent) {
    switch (parent->nchildren) {
    case 1:
        assert(child == parent->child[0]);
        parent->child[0] = NULL;
        parent->nchildren = 0;
        break;
    case 2:
        if(parent->child[1] == child)
            parent->child[1] = NULL;
        else {
            assert(parent->child[0] == child);
            parent->child[0] = parent->child[1];
            parent->child[1] = NULL;
        }
        parent->nchildren = 1;
        break;
    default:
        fprintf(stderr, "%s:%d: illegal number of children: %d\n",
                __FILE__, __LINE__, parent->nchildren);
        exit(EXIT_FAILURE);
    }
}

/// Free node and descendants.
void PopNode_free(PopNode * self) {
    if(self == NULL)
        return;

    // Free children first. Calls to PopNode_free will decrement
    // self->nchildren.
    while(self->nchildren > 0) {
        int i = self->nchildren - 1;
        PopNode_free(self->child[i]);
    }

    // unlink current node from its parents
    while(self->nparents > 0) {
        int i = self->nparents - 1;
        unlink_child(self, self->parent[i]);
        self->nparents -= 1;
    }

    free(self->label);

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

/// Duplicate a network of nodes, returning a pointer to the
/// root of the duplicate network. On entry, ppm should be an empty
/// hashmap.
PopNode *PopNode_dup(PopNode * old_root, PtrPtrMap * ppm) {
    assert(old_root);
    assert(0 == PtrPtrMap_size(ppm));
    PopNode_clear(old_root);

    // Traverse the old network, duplicating each node and
    // storing the duplicates in ppm, which maps old nodes to
    // new ones.
    PopNode_duplicate_nodes(old_root, ppm);

    // Put the old nodes into an array.
    unsigned nnodes = PtrPtrMap_size(ppm);
    void *old_nodes[nnodes];
    int status = PtrPtrMap_keys(ppm, nnodes, old_nodes);
    if(status) {
        fprintf(stderr, "%s:%d: buffer overflow\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    PopNode *node, *new_root = NULL;

    // Connect each node to its parents and children,
    // and identify the root of the duplicated network.
    for(unsigned i = 0; i < nnodes; ++i) {
        PopNode *old = old_nodes[i];
        PopNode *new = PtrPtrMap_get(ppm, old, &status);
        assert(status == 0);

        // root is the node with no parents
        if(old->nparents == 0) {
            assert(new_root == NULL);
            new_root = new;
        }
        // connect new node to its parents
        if(old->nparents > 0) {
            node = PtrPtrMap_get(ppm, old->parent[0], &status);
            assert(status == 0);
            new->parent[0] = node;
        }
        if(old->nparents > 1) {
            node = PtrPtrMap_get(ppm, old->parent[1], &status);
            assert(status == 0);
            new->parent[1] = node;
        }
        // connect new node to its children
        if(old->nchildren > 0) {
            node = PtrPtrMap_get(ppm, old->child[0], &status);
            assert(status == 0);
            new->child[0] = node;
        }
        if(old->nchildren > 1) {
            node = PtrPtrMap_get(ppm, old->child[1], &status);
            assert(status == 0);
            new->child[1] = node;
        }
    }

    return new_root;
}

/// Traverse tree, making a duplicate of each node, and putting
/// the duplicates into a hash map (called ppm) in which the old
/// node is the key and the new duplicate is the value associated
/// with that key.
static void PopNode_duplicate_nodes(PopNode * old, PtrPtrMap * ppm) {
    assert(old);
    if(old->visited)
        return;

    if(old->nsamples > 0) {
        fprintf(stderr, "%s:%d: Must call PopNode_clear before %s\n",
                __FILE__, __LINE__, __func__);
        exit(EXIT_FAILURE);
    }
    PopNode *new = memdup(old, sizeof(*old));
    CHECKMEM(new);
    old->visited = 1;
#ifndef NDEBUG
    int status = PtrPtrMap_insert(ppm, old, new);
    assert(status == 0);
#else
    (void) PtrPtrMap_insert(ppm, old, new);
#endif

    new->label = strdup(old->label);
    CHECKMEM(new->label);

    if(old->nchildren > 0)
        PopNode_duplicate_nodes(old->child[0], ppm);
    if(old->nchildren > 1)
        PopNode_duplicate_nodes(old->child[1], ppm);
}

#ifdef TEST

#  include "ptrqueue.h"
#  include "param.h"
#  include <string.h>
#  include <assert.h>
#  include <time.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

int main(int argc, char **argv) {

    int status, verbose = 0;

    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xpopnode [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, 1u);

    PtrQueue *fixedQ = PtrQueue_new();
    PtrQueue *freeQ = PtrQueue_new();
    PtrQueue *constrQ = PtrQueue_new();

    Param *par;

    par = Param_new("zero", 0.0, 0.0, 0.0, TIME | FIXED, NULL);
    PtrQueue_push(fixedQ, par);

    par = Param_new("one", 1.0, 1.0, 1.0, TWON | FIXED, NULL);
    PtrQueue_push(fixedQ, par);

    par = Param_new("Nab", 3.0, 0.0, 100.0, TWON | FREE, NULL);
    PtrQueue_push(freeQ, par);

    par = Param_new("Tab", 2.0, 0.0, 100.0, TIME | FREE, NULL);
    PtrQueue_push(freeQ, par);

    par = Param_new("Tmig", 1.0, 1.0, 1.0, TIME | FIXED, NULL);
    PtrQueue_push(fixedQ, par);

    par = Param_new("mix", 0.02, 0.02, 0.02, MIXFRAC | FIXED, NULL);
    PtrQueue_push(fixedQ, par);

    par = Param_new("Nabc", 3.0, 0.0, 100.0, TWON | FREE, NULL);
    PtrQueue_push(freeQ, par);

    par = Param_new("Tabc", 4.0, -DBL_MAX, DBL_MAX, TIME | CONSTRAINED,
                    "Tab + Nab*Nabc");
    PtrQueue_push(constrQ, par);

    ParStore *ps = ParStore_new(fixedQ, freeQ, constrQ);

    if(verbose)
        ParStore_print(ps, stderr);

    Bounds bnd = {
        .lo_twoN = 0.0,
        .hi_twoN = 1e12,
        .lo_t = 0,
        .hi_t = 1e10
    };

    PopNode *a, *b, *b2, *c, *c2, *ab, *abc;
    int ni, ti, mi;
    tipId_t ida = 0;
    tipId_t idb = 1;
    tipId_t idc = 2;
    Gene *ga = Gene_new(ida);
    Gene *gb = Gene_new(idb);
    Gene *gc = Gene_new(idc);

    ni = ParStore_getIndex(ps, "one");
    assert(ni >= 0);
    ti = ParStore_getIndex(ps, "zero");
    assert(ti >= 0);

    a = PopNode_new(ni, ti, ps, "a");
    assert(a);

    b = PopNode_new(ni, ti, ps, "b");
    assert(b);

    c = PopNode_new(ni, ti, ps, "c");
    assert(c);

    ti = ParStore_getIndex(ps, "Tmig");
    assert(ti >= 0);
    b2 = PopNode_new(ni, ti, ps, "b2");
    c2 = PopNode_new(ni, ti, ps, "c2");
    assert(b2);
    assert(c2);

    ni = ParStore_getIndex(ps, "Nab");
    assert(ni >= 0);
    ti = ParStore_getIndex(ps, "Tab");
    assert(ti >= 0);
    ab = PopNode_new(ni, ti, ps, "ab");
    assert(ab);

    ni = ParStore_getIndex(ps, "Nabc");
    assert(ni >= 0);
    ti = ParStore_getIndex(ps, "Tabc");
    assert(ti >= 0);
    abc = PopNode_new(ni, ti, ps, "abc");
    assert(abc);

    status = PopNode_addChild(ab, a);
    assert(status == 0);

    mi = ParStore_getIndex(ps, "mix");
    assert(mi >= 0);
    status = PopNode_mix(b, mi, c2, b2, ps);

    status = PopNode_addChild(c2, c);
    assert(status == 0);

    status = PopNode_addChild(ab, b2);
    assert(status == 0);

    status = PopNode_addChild(abc, ab);
    assert(status == 0);

    status = PopNode_addChild(abc, c2);
    assert(status == 0);

    if(verbose) {
        PopNode_printShallow(abc, stdout);
        PopNode_printShallow(ab, stdout);
        PopNode_printShallow(a, stdout);
        PopNode_printShallow(b, stdout);
        PopNode_printShallow(c, stdout);
    }

    assert(PopNode_isClear(abc));

    PopNode_transferSample(a, ga);
    PopNode_transferSample(b, gb);
    PopNode_transferSample(c, gc);

    assert(!PopNode_isClear(abc));

    long unsigned event_counter = 0;

    PopNode_unvisit(abc);
    Gene *root = PopNode_coalesce(abc, rng, &event_counter);
    assert(root != NULL);

    assert(!PopNode_isClear(abc));
    PopNode_clear(abc);
    assert(PopNode_isClear(abc));

    assert(abc == PopNode_root(a));
    assert(abc == PopNode_root(b));
    assert(abc == PopNode_root(b2));
    assert(abc == PopNode_root(c));
    assert(abc == PopNode_root(c2));
    assert(abc == PopNode_root(ab));
    assert(abc == PopNode_root(abc));

    PopNode_clear(abc);
    assert(PopNode_feasible(abc, bnd, verbose));

    PtrPtrMap *ppm = PtrPtrMap_new(512);
    PopNode *duproot = PopNode_dup(abc, ppm);
    CHECKMEM(duproot);
    assert(PopNode_feasible(duproot, bnd, verbose));
    Gene_free(root);
    PtrPtrMap_free(ppm);

    assert(PopNode_equals(abc, duproot));

    PopNode_free(abc);
    PopNode_free(duproot);

    unitTstResult("PopNode", "OK");

    PtrQueue_free(freeQ);
    PtrQueue_free(fixedQ);
    PtrQueue_free(constrQ);

    ParStore_free(ps);
    gsl_rng_free(rng);

    return 0;

}
#endif
