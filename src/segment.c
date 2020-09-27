/**
 * @file segment.c
 * @author Alan R. Rogers
 * @brief A single segment of a population tree.
 *
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

//#define VERBOSE

int segnum = 0;

int verbosity = 0;

#include "binary.h"
#include "branchtab.h"
#include "comb.h"
#include "error.h"
#include "idset.h"
#include "idsetset.h"
#include "matcoal.h"
#include "misc.h"
#include "parstore.h"
#include "ptrlst.h"
#include "ptrptrmap.h"
#include "segment.h"
#include "setpart.h"
#include <stdlib.h>
#include <string.h>

// Site pattern representing the union of all samples.
tipId_t union_all_samples = 0;

typedef struct MigDat MigDat;
typedef struct CombDat CombDat;
typedef struct SetPartDat SetPartDat;

// For assigning initial tipId_t values within segments.
//static tipId_t currTipId = 0;

// Data manipulated by visitComb function
struct CombDat {
    long double contrib; // Pr[site pattern]*E[len of interval]
    IdSetSet *d;
    BranchTab *branchtab;
    int dosing;    // do singleton site patterns if nonzero
};

// Data manipulated by migrate function
struct MigDat {
    int nMigrants, nNatives;
    long double mig_pr; // probability of this migration outcome
    PtrLst *migrants, *natives, *a;
    unsigned mig_event, mig_outcome;
};

// Data manipulated by visitSetPart function.
struct SetPartDat {
    unsigned nparts; // Number of parts in partition

    long double prior; // prob of k ancestors given n descendants
    long double lnconst;  // log of constant in Durrett's theorem 1.5
    long double elen;     // E[len of interval]

    // a is a list of IdSet objects for sets of ancestors
    PtrLst *a;

    // descendants
    IdSetSet *d;

    BranchTab *branchtab;
    int dosing;      // do singleton site patterns if nonzero
};

// One segment of a population network. This version works
// with MCTree.
struct Segment {
    int            segnum; // for debugging
    int            nparents, nchildren, nsamples;
    int            visited;     // for traversal algorithm
    double    twoN;        // ptr to current pop size
    double    start, end;  // duration of this Segment
    double    mix;         // ptr to frac of pop derived from parent[1]

    // indices into ParStore array
    int twoN_i, start_i, end_i, mix_i;

    struct Segment *parent[2];
    struct Segment *child[2];

    tipId_t sample[MAXSAMP]; // array of length nsamples

    // Waiting rooms.  Each child loads IdSet objects into one of two
    // waiting rooms. w[i] is the i'th waiting room, where i is in
    // [0,1]. w[i] is the waiting room for child i. The number of
    // waiting rooms equals the number of children. w[i][j] is a list
    // of IdSet objects, each of which contains j tipId_t values.
    // wdim[i] is the dimension of the i'th waiting room, so 0 <= j <
    // wdim[i]. If there are no waiting rooms (i.e. if nchildren=0),
    // then the descendants at the beginning of the segment consist
    // only of those in array "sample". If there is one waiting room,
    // then the ids in "sample" are added to each of the sets in w[0],
    // and the resulting list of IdSet objects represents the state at
    // the beginning of the segment. If there are two waiting rooms,
    // then we begin with all pairs of IdSet objects in w[0] and
    // w[1]. To each pair, we add the ids in "samples" to obtain the
    // list of IdSet objects at the beginning of the segment.
    int wdim[2];
    IdSetSet *w[2][MAXSAMP+1];

    // Dimension of array d. d[i] holds IdSet objects with i lineages.
    // Thus d[0] refers to empty sets. The largest possible set has
    // dim-1 lineages.
    int dim;

    // d[i] is a pointer to an IdSetSet, which is allocated within
    // Segment_coalesce. This holds IdSet objects that contain i
    // tipId_t values, each representing a descendant at the recent
    // end of the Segment. The dimension of array d is self->dim.
    IdSetSet **d;
};

static IdSetSet **get_descendants1(int dim, IdSetSet **w, int nsamples,
                                   tipId_t *sample, int *newdim);
static IdSetSet **get_descendants2(int dim0, IdSetSet **w0,
                                   int dim1, IdSetSet **w1,
                                   int nsamples, tipId_t *sample, int *newdim);
static void mv_idsets_to_parent(Segment *self, int ipar, PtrLst **a);
static void coalescent_interval_length(int n, long double elen[n-1],
                                       long double eig[n],
                                       long double v);
static void project(int n, long double pr[n], long double eig[n-1]);
int    visitComb(int d, int ndx[d], void *data);
int    visitSetPart(unsigned n, unsigned a[n], void *data);
int    visitMig(int nmig, int *migndx, void *data);
static void  unlink_child(Segment *child, Segment *parent);
static int Segment_coalesceFinite(Segment *self, int dosing,
                                  BranchTab *branchtab);
static int Segment_coalesceInfinite(Segment *self, long double v, int dosing,
                                    BranchTab *branchtab);
static void Segment_duplicate_nodes(Segment *old, PtrPtrMap *ppm);
static int Segment_equals_r(Segment *a, Segment *b);
static int self_ndx(Segment *self, Segment *parent);
static int w_isempty(int dim, IdSetSet **w);
static void migrate(PtrLst *migrants, PtrLst *natives, PtrLst *sets,
                    int nmig, int *migndx, int nnat, int *natndx,
                    MigDat *md);
static void mv_to_waiting_room(Segment *self, PtrLst *src, int ipar,
                               int nlin);
void Segment_print_d(Segment *self, const char *func, int line);
static int tipidcmp(const void *vx, const void *vy);

// Return index of self among children of parent
static int self_ndx(Segment *self, Segment *parent) {
    if(self == parent->child[0])
        return 0;
    else if(self == parent->child[1])
        return 1;
    else {
        fprintf(stderr,"%s:%d: self not child of parent\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
}

/// Remove child from parent
static void unlink_child(Segment *child, Segment *parent) {
    switch(parent->nchildren) {
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
        fprintf(stderr,"%s:%d: illegal number of children: %d\n",
                __FILE__,__LINE__, parent->nchildren);
        exit(EXIT_FAILURE);
    }
}

void *Segment_new(int twoN_i, int start_i, ParStore *ps) {
    
    Segment *self = malloc(sizeof(*self));
    CHECKMEM(self);

    memset(self, 0, sizeof(*self));

    self->segnum = segnum++;  // debug

    self->twoN_i = twoN_i;
    self->start_i = start_i;
    self->end_i = -1;
    self->mix_i = -1;

    self->twoN = ParStore_getVal(ps, twoN_i);
    self->start = ParStore_getVal(ps, start_i);
    self->end = INFINITY;
    self->mix = 0.0;

    for(int i=0; i <= MAXSAMP; ++i) {
        self->w[0][i] = IdSetSet_new(0);
        self->w[1][i] = IdSetSet_new(0);
    }

    // other fields are not initialized by Segment_new.

    return self;
}

void Segment_free(Segment *self) {
    if(self == NULL)
        return;

    // Free children first. Calls to Segment_free will decrement
    // self->nchildren.
    while(self->nchildren) {
        int i = self->nchildren - 1;
        Segment_free(self->child[i]);
    }

    // unlink current node from its parents
    while(self->nparents > 0) {
        int i = self->nparents - 1;
        unlink_child(self, self->parent[i]);
        self->nparents -= 1;
    }

    for(int i=0; i <= MAXSAMP; ++i) {
        IdSetSet_free_deep(self->w[0][i]);
        IdSetSet_free_deep(self->w[1][i]);
    }

    if(self->d)
        free(self->d);

    free(self);
}

/// Traverse network, removing segments with no children and
/// no samples. These do not contribute to the calculation.
void Segment_prune(Segment *self) {
    if(self->nchildren > 1)
        Segment_prune(self->child[1]);
    if(self->nchildren > 0)
        Segment_prune(self->child[0]);

    if(self->nchildren == 0 && self->nsamples == 0) {
        if(self->nparents == 2)
            unlink_child(self, self->parent[1]);
        if(self->nparents > 0)
            unlink_child(self, self->parent[0]);
        Segment_free(self);
    }
}

/// Set all "visited" flags to false.
void Segment_unvisit(Segment *self) {
    if(self->nchildren > 0)
        Segment_unvisit(self->child[0]);
    if(self->nchildren > 1)
        Segment_unvisit(self->child[1]);
    self->visited = 0;
}

/// Add a new sample to a Segment. This is called once for each
/// of the samples specified in the .lgo file.
void Segment_newSample(Segment * self, unsigned ndx) {
    assert(1 + self->nsamples < MAXSAMP);
    assert(ndx < 8 * sizeof(tipId_t));

    static const tipId_t one = 1;
    self->sample[self->nsamples] = one << ndx;
    self->nsamples += 1;
}

int Segment_addChild(void * vparent, void * vchild) {
    Segment *parent = vparent;
    Segment *child = vchild;
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
                    __FILE__, __func__, __LINE__,
                    child->end_i, parent->start_i);
            fprintf(stderr, "  child->end=%lg != %lg = parent->start\n",
                    child->end, parent->start);
            return DATE_MISMATCH;
    }

    parent->child[parent->nchildren] = child;
    child->parent[child->nparents] = parent;
    ++parent->nchildren;
    ++child->nparents;
    Segment_sanityCheck(parent, __FILE__, __LINE__);
    Segment_sanityCheck(child, __FILE__, __LINE__);
    return 0;
}

/// Check sanity of Segment
void Segment_sanityCheck(Segment * self, const char *file, int lineno) {
#ifndef NDEBUG
    REQUIRE(self != NULL, file, lineno);
    REQUIRE(self->twoN > 0.0, file, lineno);
    REQUIRE(self->start >= 0.0, file, lineno);
    if(self->end_i >= 0)
        REQUIRE(self->nparents > 0, file, lineno);
    REQUIRE(self->start <= self->end, file, lineno);
    switch(self->nparents) {
    case 0:
        REQUIRE(self->end_i == -1, file, lineno);
        REQUIRE(self->mix_i == -1, file, lineno);
        break;
    case 1:
        REQUIRE(self->end_i >= 0, file, lineno);
        REQUIRE(self->mix_i == -1, file, lineno);
        break;
    case 2:
        REQUIRE(self->end_i >= 0, file, lineno);
        REQUIRE(self->mix_i >= 0, file, lineno);
        break;
    default:
        fprintf(stderr,"%s:%d: bad number of parents: %d\n",
                file, lineno, self->nparents);
        exit(EXIT_FAILURE);
    }
    REQUIRE(self->nsamples >= 0, file, lineno);
    REQUIRE(self->nsamples < MAXSAMP, file, lineno);
    REQUIRE(self->nchildren >= 0, file, lineno);
    REQUIRE(self->nchildren <= 2, file, lineno);
    if(self->nchildren > 0)
        Segment_sanityCheck(self->child[0], file, lineno);
    if(self->nchildren > 1)
        Segment_sanityCheck(self->child[1], file, lineno);

    REQUIRE(no_shared_bits(self->nsamples, self->sample), file, lineno);
    REQUIRE(self->wdim[0] <= MAXSAMP+1, file, lineno);
    REQUIRE(self->wdim[1] <= MAXSAMP+1, file, lineno);
    IdSet *id;
    if(self->nchildren > 0) {
        for(int i=0; i < self->wdim[0]; ++i) {
            IdSetSet_rewind(self->w[0][i]);
            for(id=IdSetSet_next(self->w[0][i]);
                id;
                id=IdSetSet_next(self->w[0][i])) {
                
                IdSet_sanityCheck(id, __FILE__, __LINE__);
                
            }
        }
    }
    if(self->nchildren > 1) {
        for(int i=0; i < self->wdim[1]; ++i) {
            IdSetSet_rewind(self->w[1][i]);
            for(id=IdSetSet_next(self->w[1][i]);
                id;
                id=IdSetSet_next(self->w[1][i])) {
                
                IdSet_sanityCheck(id, file, lineno);
                
            }
        }
    }
    if(self->d != NULL) {
        for(int i=0; i < self->dim; ++i) {
            IdSetSet_rewind(self->d[i]);
            for(id = IdSetSet_next(self->d[i]);
                id;
                id = IdSetSet_next(self->d[i])) {

                IdSet_sanityCheck(id, file, lineno);

            }
        }
    }
#endif
}

int Segment_mix(void * vchild, int mix_i, void * vintrogressor,
                void * vnative, ParStore *ps) {
    Segment *child = vchild, *introgressor = vintrogressor,
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
            fprintf(stderr,"%s:%s:%d: Date mismatch\n"
                    "  child->end_i=%d != %d=introgressor->start_i\n",
                    __FILE__, __func__, __LINE__,
                    child->end_i, introgressor->start_i);
            fprintf(stderr,"  child->end=%lg != %lg=introgressor->start\n",
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
    Segment_sanityCheck(child, __FILE__, __LINE__);
    Segment_sanityCheck(introgressor, __FILE__, __LINE__);
    Segment_sanityCheck(native, __FILE__, __LINE__);
    return 0;
}

/// Find root of population network, starting from given node.
void    *Segment_root(void * vself) {
    Segment *self = vself, *r0, *r1;
    assert(self);
    switch (self->nparents) {
    case 0:
        return self;
        break;
    case 1:
        return Segment_root(self->parent[0]);
        break;
    case 2:
        r0 = Segment_root(self->parent[0]);
        r1 = Segment_root(self->parent[1]);
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

/// Return 1 if parameters satisfy inequality constraints, or 0 otherwise.
int Segment_feasible(const Segment * self, Bounds bnd, int verbose) {
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
        if(0 == Segment_feasible(self->child[i], bnd, verbose))
            return 0;
    }
    return 1;
}

/// Duplicate a network of nodes, returning a pointer to the
/// root of the duplicate network. On entry, ppm should be an empty
/// hashmap.
Segment *Segment_dup(Segment *old_root, PtrPtrMap *ppm) {
    assert(old_root);
    assert(0 == PtrPtrMap_size(ppm));

    // Traverse the old network, duplicating each node and
    // storing the duplicates in ppm, which maps old nodes to
    // new ones.
    Segment_unvisit(old_root);
    Segment_duplicate_nodes(old_root, ppm);

    // Put the old nodes into an array.
    unsigned nnodes = PtrPtrMap_size(ppm);
    void *old_nodes[nnodes];
    int status = PtrPtrMap_keys(ppm, nnodes, old_nodes);
    if(status) {
        fprintf(stderr,"%s:%d: buffer overflow\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    Segment *node, *new_root=NULL;

    // Connect each node to its parents and children,
    // identify the root of the duplicated network,
    // initialize nu->d and nu->w.
    for(unsigned i=0; i < nnodes; ++i) {
        Segment *old = old_nodes[i];
        Segment *nu = PtrPtrMap_get(ppm, old, &status);
        assert(status == 0);

        // root is the node with no parents
        if(old->nparents == 0) {
            assert(new_root == NULL);
            new_root = nu;
        }

        // connect new node to its parents
        if(old->nparents > 0) {
            node = PtrPtrMap_get(ppm, old->parent[0], &status);
            assert(status == 0);
            nu->parent[0] = node;
        }
        if(old->nparents > 1) {
            node = PtrPtrMap_get(ppm, old->parent[1], &status);
            assert(status == 0);
            nu->parent[1] = node;
        }

        // connect new node to its children
        if(old->nchildren > 0) {
            node = PtrPtrMap_get(ppm, old->child[0], &status);
            assert(status == 0);
            nu->child[0] = node;
        }
        if(old->nchildren > 1) {
            node = PtrPtrMap_get(ppm, old->child[1], &status);
            assert(status == 0);
            nu->child[1] = node;
        }

        nu->d = NULL;
        nu->wdim[0] = old->wdim[0];
        nu->wdim[1] = old->wdim[1];

        for(int j=0; j <= MAXSAMP; ++j) {
            nu->w[0][j] = IdSetSet_new(0);
            nu->w[1][j] = IdSetSet_new(0);
        }
    }

    return new_root;
}

void     Segment_print(FILE * fp, void * vself, int indent) {
    Segment *self = vself;
    for(int i = 0; i < indent; ++i)
        fputs("   ", fp);
    fprintf(fp, "%p twoN=%lf t=(%lf,", self, self->twoN, self->start);
    fprintf(fp, "%lf)", self->end);
    if(self->nchildren)
        fprintf(fp, " nchildren=%d", self->nchildren);
    if(self->nsamples)
        fprintf(fp, " nsamples=%d", self->nsamples);
    putc('\n', fp);

    for(int i = 0; i < self->nchildren; ++i)
        Segment_print(fp, self->child[i], indent + 1);
}

/// Visit a combination
int visitComb(int d, int ndx[d], void *data) {
    assert(d>0);
    assert(union_all_samples != 0);
    CombDat *dat = (CombDat *) data;
    IdSet *ids;

    IdSetSet_rewind(dat->d);
    while( (ids = IdSetSet_next(dat->d)) != NULL) {
        
        // sitepat is the union of the current set of descendants, as
        // described in ndx.
        tipId_t sitepat = 0;
        for(int j=0; j < d; ++j) {
            assert(ndx[j] < ids->nIds);
            assert(ndx[j] >= 0);
            assert(ids->tid[ndx[j]]);
            sitepat |= ids->tid[ndx[j]];
        }
        
        assert(sitepat > 0);

        // Skip singletons unless data->dosing is nonzero
        if(!dat->dosing && isPow2(sitepat))
            continue;

        if(sitepat == union_all_samples)
            continue;

#ifdef VERBOSE        
        fprintf(stderr,"%s:%d: adding %Lg to pattern o%o\n",
                __FILE__,__LINE__, IdSet_prob(ids) * dat->contrib,
                sitepat);
#endif        
      
        // Increment BranchTab entry for current sitepat value.
        BranchTab_add(dat->branchtab, sitepat,
                      IdSet_prob(ids) * dat->contrib);
    }
    return 0;
}

/// Visit a set partition. n is the number of descendants, a[i] is the
/// index of the ancestor of the i'th descendant.
int visitSetPart(unsigned n, unsigned a[n], void *data) {
    assert(union_all_samples != 0);
    SetPartDat *vdat = data;
    int status=0;

    // Calculate the size, c[i], of each part. In other words,
    // the number of descendants of each ancestor.
    unsigned k=vdat->nparts, c[k];
    memset(c, 0, k*sizeof(c[0]));
    for(int i=0; i<n; ++i) {
        assert(a[i] < k);
        ++c[a[i]];
    }

    tipId_t sitepat[k];

    long double p = probPartition(k, c, vdat->lnconst);

#ifdef VERBOSE    
    fprintf(stderr,"%s:%s:%d: set prob=%Lg prior=%Lg\n",
            __FILE__,__func__,__LINE__, p, vdat->prior);
#endif    


    // Loop over IdSet objects
    // a set of n descendants.
    IdSetSet_rewind(vdat->d);
    for(IdSet *descendants = IdSetSet_next(vdat->d);
        descendants;
        descendants = IdSetSet_next(vdat->d)) {

        memset(sitepat, 0, k*sizeof(tipId_t));

        assert(n == IdSet_nIds(descendants));

        // Loop over descendants, creating a sitepat for each
        // ancestor. a[i] is the index of the ancestor of the i'th
        // descendant. sitepat[j] is the site pattern of the j'th
        // ancestor.
        for(int j=0; j<n; ++j)
            sitepat[a[j]] |= descendants->tid[j];

        // Loop over ancestors, i.e. over site patterns, adding
        // to the corresponding entry in BranchTab.
        for(int j=0; j<k; ++j) {

            if(sitepat[j] == union_all_samples)
                continue;

#ifdef VERBOSE            
            fprintf(stderr,"%s:%d: adding %Lg=%Lg*%Lg*%Lg to pattern o%o\n",
                    __FILE__,__LINE__,
                    p * vdat->elen * IdSet_prob(descendants),
                    p, vdat->elen, IdSet_prob(descendants),
                    sitepat[j]);
            IdSet_print(descendants, stderr);
#endif            
      
            BranchTab_add(vdat->branchtab, sitepat[j],
                          p * vdat->elen * IdSet_prob(descendants));
        }

        // Add the current set partition to the list of ancestral
        // states.

        qsort(&sitepat[0], (size_t) k, sizeof(sitepat[0]), tipidcmp);

#ifdef VERBOSE        
        fprintf(stderr,"%s:%s:%d: new IdSet pr=%Lg * %Lg\n",
                __FILE__,__func__,__LINE__, p, descendants->p);
#endif        
        IdSet *ancestors = IdSet_new(k, sitepat,
                                     p * vdat->prior * descendants->p);
        IdSet_copyMigOutcome(ancestors, descendants);

#ifndef NDEBUG        
        IdSet_sanityCheck(ancestors, __FILE__, __LINE__);
#endif        
        status = PtrLst_push(vdat->a, ancestors);
        if(status) {
            fprintf(stderr,"%s:%d can't push ancestors; status=%d\n",
                    __FILE__,__LINE__, status);
            exit(EXIT_FAILURE);
        }
    }

    return status;
}

/// Visit a combination defining migrants. migndx has the 0-based
/// indices of the current set of migrants. Its length is nmig, which
/// may be zero.
int visitMig(int nmig, int *migndx, void *data) {
    MigDat *mdat = (MigDat *) data;
    int nnat = mdat->nNatives;
    int i;

    assert(nmig == mdat->nMigrants);

    // Array natndx has size nnat. I allocate an extra entry in case
    // nnat == 0.
    int natndx[1 + nnat];
    
    // set natndx equal to complement of migndx
    int next = 0, j = 0;
    for(i=0; i<nmig; ++i) {
        while(next < migndx[i])
            natndx[j++] = next++;
        next = migndx[i] + 1;
    }
    while(j < nnat)
        natndx[j++] = next++;

#if 1
    fprintf(stderr,"%s:%d: nmig=%d nnat=%d event_outcome=%d_%d: ",
            __func__,__LINE__, nmig, nnat, mdat->mig_event,
            mdat->mig_outcome);
    for(int ii=0; ii<nmig; ++ii)
        fprintf(stderr," %d", migndx[ii]);
    fputs(" :", stderr);
    for(int ii=0; ii<nnat; ++ii)
        fprintf(stderr," %d", natndx[ii]);
    putc('\n', stderr);
#endif

    migrate(mdat->migrants, mdat->natives, mdat->a, nmig, migndx,
            nnat, natndx, mdat);

    return 0;
}

/**
 * On entry, "sets" is a vector of IdSet objects, which will be copied,
 * and the copies allocated among migrants and natives. Each copy will
 * acquire a MigOutcome label, which is constructed using mig_event,
 * mig_outcome, and mig_pr. Each set within "sets" should have the
 * same length (number of tipId_t values), and that length should
 * correspond to the indices in arrays "migndx" and "natndx". The
 * first of these arrays holds the indices of the migrant tipId_t
 * values within each set, and the second holds indices of natives.
 * On return, "migrants" and "natives" contain IdSet objects of
 * migrants and natives.
 */
static void migrate(PtrLst *migrants, PtrLst *natives, PtrLst *sets,
                    int nmig, int *migndx, int nnat, int *natndx,
                    MigDat *mdat) {

    IdSet *set;
    
    // Extra entry guards against problems of zero length.
    tipId_t migid[1 + nmig], natid[1 + nnat];

    // number of sets of ancestors
    int nSets = PtrLst_length(sets);

    fprintf(stderr,"%s:%d: nSets=%d\n", __func__,__LINE__, nSets);

    PtrLst_rewind(sets);
    while( (set = PtrLst_next(sets)) != NULL ) {
        
        fprintf(stderr,"%s:%d:", __func__,__LINE__);
        for(int i=0; i < nmig+nnat; ++i)
            fprintf(stderr," %o", set->tid[i]);
        fputs(" : ", stderr);
        MigOutcome_print(set->mig, stderr);
        putc('\n', stderr);

        IdSet_sanityCheck(set, __FILE__, __LINE__);

        assert(IdSet_nIds(set) == nmig + nnat);

        for(int i=0; i<nmig; ++i)
            migid[i] = set->tid[migndx[i]];

#if 1        
        fprintf(stderr,"%s:%d: %d_%d: migrants:",
                __func__,__LINE__, mdat->mig_event, mdat->mig_outcome);
        for(int i=0; i<nmig; ++i)
            fprintf(stderr," %o", migid[i]);
        putc('\n', stderr);
#endif        

        for(int i=0; i < nnat; ++i)
            natid[i] = set->tid[natndx[i]];

#if 1        
        fprintf(stderr,"%s:%d: %d_%d: natives:",
                __func__,__LINE__, mdat->mig_event, mdat->mig_outcome);
        for(int i=0; i<nnat; ++i)
            fprintf(stderr," %o", natid[i]);
        putc('\n', stderr);
#endif

        // Create IdSet objects for migrants and natives
        IdSet *mig = IdSet_new(nmig, migid, set->p);
        IdSet_addMigEvent(mig, mdat->mig_event, mdat->mig_outcome,
                          mdat->mig_pr);
#ifndef NDEBUG        
        IdSet_sanityCheck(mig, __FILE__, __LINE__);
#endif        
        PtrLst_push(migrants, mig);

        IdSet *nat = IdSet_new(nnat, natid, set->p);
        IdSet_addMigEvent(nat, mdat->mig_event, mdat->mig_outcome,
                          mdat->mig_pr);
        ++mdat->mig_outcome;
#ifndef NDEBUG        
        IdSet_sanityCheck(nat, __FILE__, __LINE__);
#endif        
        PtrLst_push(natives, nat);
    }
}

/// Traverse tree, making a duplicate of each node, and putting
/// the duplicates into a hash map (called ppm) in which the old
/// node is the key and the new duplicate is the value associated
/// with that key.
static void Segment_duplicate_nodes(Segment *old, PtrPtrMap *ppm) {
    assert(old);
    if(old->visited)
        return;

    Segment *new = memdup(old, sizeof(*old));
    CHECKMEM(new);
    old->visited = 1;
    int status = PtrPtrMap_insert(ppm, old, new);
    assert(status==0);
    if(old->nchildren > 0)
        Segment_duplicate_nodes(old->child[0], ppm);
    if(old->nchildren > 1)
        Segment_duplicate_nodes(old->child[1], ppm);

    for(int i=0; i <= MAXSAMP; ++i) {
        new->w[0][i] = IdSetSet_new(0);
        new->w[1][i] = IdSetSet_new(0);
    }
}

int Segment_equals(Segment *a, Segment *b) {
    Segment_unvisit(a);
    Segment_unvisit(b);
    return Segment_equals_r(a, b);
}

static int Segment_equals_r(Segment *a, Segment *b) {
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
    for(int i=0; i < a->nchildren; ++i) {
        if(!Segment_equals(a->child[i], b->child[i]))
            return 0;
    }
    return 1;
}

static int Segment_coalesceFinite(Segment *self, int dosing,
                                  BranchTab *branchtab) {

    assert(union_all_samples != 0);
    Segment_print_d(self, __func__,__LINE__);

    /*
      pr[i] = prob of i+1 lineages
      elen[i] = expected length of interval w/ i+1 lineages
      a[i] = list of ancestral IdSet objects w/ i lineages
     */
    assert(self->dim > 0);
    long double pr[self->dim], elen[self->dim];
    long double v = (self->end - self->start) / self->twoN;
    int n, i, k, iself, status=0;

    assert(isfinite(v));
    assert(v >= 0.0);

    // Array of lists of ancestors. a[k] is the list for sets of k
    // ancestors. The case of k=0 arises from migration events in
    // which all migrants went the other way. When these empty IdSet
    // objects combine with full objects from a sister branch, they
    // contribute information about migration history.
    PtrLst *a[self->dim];
    for(i=0; i < self->dim; ++i)
        a[i] = PtrLst_new();

    SetPartDat sd = {.branchtab = branchtab,
                     .dosing = dosing,
    };

    // Cases of 0 or 1 lineages
    for(n=0; n <= 1 && n < self->dim; ++n) {
        int nIds = IdSetSet_size(self->d[n]);

        if(0 == nIds)
            continue;

        IdSetSet_rewind(self->d[n]);
        for(IdSet *ids = IdSetSet_next(self->d[n]);
            ids;
            ids = IdSetSet_next(self->d[n])) {

            if(n==1 && ids->tid[0] != union_all_samples) {
                // Single lineage in finite Segment contributes
                // to branchtab.
                assert(IdSet_nIds(ids) == 1);

#ifdef VERBOSE                
                fprintf(stderr,"%s:%d: adding %Lg to pattern o%o\n",
                        __FILE__,__LINE__,
                        IdSet_prob(ids) * v * self->twoN,
                        ids->tid[0]);
#endif                
      
                BranchTab_add(branchtab, ids->tid[0],
                              IdSet_prob(ids) * v * self->twoN);
            }

            IdSet *new = IdSet_dup(ids);
            IdSet_sanityCheck(new, __FILE__,__LINE__);
            PtrLst_push(a[n], new);
        }
    }

    // Cases of 2..(dim-1) lineages.  Outer loop over numbers of
    // descendants.  Calculate probabilities and expected values, p[1]
    // and elen.
    for(n=2; n < self->dim; ++n) {

        // eigenvalues of transient states. eig[i] refers to
        // interval with i+2 lineages.
        long double eig[n-1];
        MatCoal_eigenvals(n-1, eig, v);

        // Calculate the expected length, elen[i], of the subinterval
        // containing i+1 lineages.
        coalescent_interval_length(n, elen, eig, v);

        // Convert from coalescent time to generations
        for(int ii=0; ii < n; ++ii)
            elen[ii] *= self->twoN;
        
        // Calculate pr[i], the probability of i+1 lineages at the
        // ancient end of the segment.
        project(n, pr, eig);

#ifdef VERBOSE        
        fprintf(stderr,"%s:%d: MatCoal pr:",__FILE__,__LINE__);
        for(int ii=0; ii<n; ++ii)
            fprintf(stderr," %d:%Lg", ii+1, pr[ii]);
        putc('\n', stderr);
        fprintf(stderr,"%s:%d: MatCoal elen:",__FILE__,__LINE__);
        for(int ii=0; ii<n; ++ii)
            fprintf(stderr," %d:%Lg", ii+1, elen[ii]);
        putc('\n', stderr);
#endif        

        sd.d = self->d[n];

        // Loop over number, k, of ancestors.  Include k=1, because
        // this is a finite segment.
        for(k=1; k <= n; ++k) {
            sd.a = a[k];
            sd.nparts = k;
            sd.prior = pr[k-1];
            sd.elen = elen[k-1];
            sd.lnconst = lnCoalConst(n, k);
            status = traverseSetPartitions(n, k, visitSetPart, &sd);
            if(status)
                return status;
        }
    }

    // Transfer IdSet objects to parental waiting rooms.
    if(self->nparents == 1) {
        if(self->segnum == 4) {
            fprintf(stderr,"%s:%d: verbosity on for segnum %d\n",
                    __func__,__LINE__, self->segnum);
            verbosity = 1;
        }
        mv_idsets_to_parent(self, 0, a);
        verbosity = 0;
    }else if(self->mix == 0.0) {

        // everyone is a native
        mv_idsets_to_parent(self, 0, a);
        
    }else if(self->mix == 1.0) {

        // everyone is a migrant
        mv_idsets_to_parent(self, 1, a);
        
    }else{
        assert(self->nparents == 2);
        assert(self->mix > 0.0);
        assert(self->mix < 1.0);

        // natives go to parent[0], migrants to parent[1]

        // Set wdim in the two parents
        for(int ipar = 0; ipar < 2; ++ipar) {
            iself = self_ndx(self, self->parent[ipar]);
            self->parent[ipar]->wdim[iself] = self->dim;
        }

        MigDat msd = {
                      .migrants = PtrLst_new(),
                      .natives = PtrLst_new(),
                      .a = NULL,
                      .mig_event = nextMigrationEvent(),
                      .mig_outcome = 0
        };
        fprintf(stderr,"%s:%d: mig_event %d is in segnum %d\n",
                __func__,__LINE__, msd.mig_event, self->segnum);

        // Loop over the number, k, of ancestors. This includes 0
        // because migration can produce empty sets (because everyone
        // migrated the other way).
        for(k=0; k < self->dim; ++k) {

            msd.a = a[k];
            if(PtrLst_length(msd.a) == 0)
                continue;

            for(int x=0; x <= k; ++x) {

                // Prob that x of k lineages are migrants is
                // binomial with index k and parameter "mix".
                long double migprob = logl(binom(k, x))
                    + x*logl(self->mix)
                    + (k-x)*logl(1.0-self->mix);
                migprob = expl(migprob);
                assert(isfinite(migprob));

                msd.mig_pr = migprob;
                msd.nMigrants = x;
                msd.nNatives = k - x;

                status = traverseComb(k, x, visitMig, &msd);
                if(status)
                    return status;

                // transfer natives
                mv_to_waiting_room(self, msd.natives, 0, k-x);

                // transfer migrants
                mv_to_waiting_room(self, msd.migrants, 1, x);
            }
            {
                IdSet *s;
                while( (s = PtrLst_pop(a[k])) != NULL)
                    IdSet_free(s);
            }
        }
        PtrLst_free(msd.migrants);
        PtrLst_free(msd.natives);
    }

    for(i=0; i < self->dim; ++i)
        PtrLst_free(a[i]);

    return status;
}

static int Segment_coalesceInfinite(Segment *self, long double v,
                                    int dosing,
                                    BranchTab *branchtab) {
    Segment_print_d(self, __func__,__LINE__);

    assert(self->dim > 0);
    long double elen[self->dim];
    int n, status=0;

    CombDat cd = {.branchtab = branchtab,
                  .dosing = dosing
    };

    // Outer loop over numbers of descendants.
    // Calculate expected branch lengths, elen.
    for(n=2; n < self->dim; ++n) {

        // In an infinite interval, eigenvalues of all transient
        // states are zero.
        long double eig[n-1];
        memset(eig, 0, (n-1)*sizeof(eig[0]));
        
        // Calculate the expected length, elen[i], of the subinterval
        // containing i+1 lineages.
        coalescent_interval_length(n, elen, eig, v);

        // Convert from coalescent time to generations
        for(int ii=0; ii < n; ++ii)
            elen[ii] *= self->twoN;
        
        if(IdSetSet_size(self->d[n]) == 0)
            continue;

        cd.d = self->d[n];

        // Loop over number, k, of ancestors.  Exclude k=1, because
        // this is an infinite segment.
        for(int k=2; k <= n; ++k) {
            
            // portion of log Qdk that doesn't involve d
            long double lnconst = logl(k) - logl(binom(n-1, k-1));

            // Within each interval, there can be ancestors
            // with 1 descendant, 2, 3, ..., n-k+1.
            for(int d=1; d <= n-k+1; ++d) {

                long double lnprob = lnconst
                    + logl(binom(n-d-1, k-2))
                    - logl(binom(n,d));

                // probability of site pattern
                cd.contrib = expl(lnprob);

                // times expected length of interval
                cd.contrib *= elen[k-1];
                
                status = traverseComb(n, d, visitComb, &cd);
                if(status)
                    return status;
            }
        }
        cd.d = NULL;
    }
    return status;
}

int Segment_coalesce(Segment *self, int dosing, BranchTab *branchtab) {
#ifdef VERBOSE    
    fprintf(stderr,"%s:%d: SEGNUM %d\n",
            __func__,__LINE__, self->segnum);
#endif    

    int status=0;

    if(self->visited)
        return 0;
    self->visited = 1;

    switch(self->nchildren) {
    case 0:
        fprintf(stderr,"%s:%d: SEGNUM %d is a tip\n",
                __func__,__LINE__, self->segnum);
        break;
    case 1:
        fprintf(stderr,"%s:%d: SEGNUM %d <- %d\n",
                __func__,__LINE__, self->segnum,
                self->child[0]->segnum);
        break;
    default:
        fprintf(stderr,"%s:%d: SEGNUM %d <- %d + %d\n",
                __func__,__LINE__, self->segnum,
                self->child[0]->segnum, self->child[1]->segnum);
        break;
    }

    if(self->nchildren > 0) {
        status = Segment_coalesce(self->child[0], dosing, branchtab);
        if(status)
            return status;
    }
    if(self->nchildren == 2 ) {
        status = Segment_coalesce(self->child[1], dosing, branchtab);
        if(status)
            return status;
    }

    // Set self->d, the array of descendants, and self->dim, the
    // dimension of this array.
    switch(self->nchildren) {
    case 0:
        self->d = get_descendants1(0, NULL,
                                   self->nsamples, self->sample,
                                   &self->dim);
        if(self->dim == 0)
            return 0;
        break;
    case 1:
        self->d = get_descendants1(self->wdim[0], self->w[0],
                                   self->nsamples, self->sample,
                                   &self->dim);
        if(self->dim == 0)
            return 0;
        assert(w_isempty(self->wdim[0], self->w[0]));
        break;
    default:
        assert(self->nchildren == 2);
        fprintf(stderr,"%s:%d: SEGNUM %d calling get_descendants2\n",
                __func__,__LINE__,self->segnum);
        self->d = get_descendants2(self->wdim[0], self->w[0],
                                   self->wdim[1], self->w[1],
                                   self->nsamples, self->sample,
                                   &self->dim);
        if(self->dim == 0)
            return 0;
        assert(w_isempty(self->wdim[0], self->w[0]));
        assert(w_isempty(self->wdim[1], self->w[1]));
        break;
    }
    assert(self->dim > 0);

#ifdef VERBOSE    
    fprintf(stderr,"%s:%s:%d: nchild=%d npar=%d d:\n",
            __FILE__,__func__,__LINE__,
            self->nchildren, self->nparents);
    for(int i=0; i < self->dim; ++i) {
        fprintf(stderr,"%2d:", i);
        IdSetSet_rewind(self->d[i]);
        IdSet *ids;
        while( (ids = IdSetSet_next(self->d[i])) != NULL ) {
            IdSet_print(ids, stderr);
        }
        putc('\n', stderr);
    }
#endif    

    if(self->end_i == -1) {
        status = Segment_coalesceInfinite(self, INFINITY, dosing, branchtab);
    }else{
        status = Segment_coalesceFinite(self, dosing, branchtab);
    }

    // Free IdSet objects of descendants
    for(int i=0; i < self->dim; ++i) {
        IdSetSet_free_deep(self->d[i]);
    }
    free(self->d);
    self->d = NULL;

    return status;
}

/// Return a newly-allocated array of IdSetSet objects. The i'th
/// object contains all the IdSet entries in w[i], with each entry
/// augmented by the entries (if any) in vector "sample". If there are
/// no entries in any of the w[i], the returned array contains a
/// single non-empty entry, with the entries of "samples". If "w" and
/// "samples" are both empty, *newdim will equal 0 and the function
/// returns NULL. Otherwise, *newdim is the dimension of the returned
/// array.
///
/// On entry, wdim is the dimension of array w, and nsamples is
/// the dimension of array sample. On return *newdim is the dimension
/// of the newly-allocated array returned by the function.
static IdSetSet **get_descendants1(int wdim, IdSetSet **w, int nsamples,
                                 tipId_t *sample, int *newdim) {
    int i, n, m, status;

    if(w) {
        while(wdim > 0 && IdSetSet_size(w[wdim-1]) == 0)
            wdim -= 1;
    }else
        wdim = 0;

    // Dimension is 1 larger than nsamples plus the largest
    // entry of w. If nsamples==wdim==0, then n==0.
    n = nsamples + wdim;
    if(nsamples > 0 && wdim == 0)
        n += 1;
    
    *newdim = n;

    if(n == 0)
        return NULL;

    IdSetSet **d = malloc(n * sizeof(IdSetSet *));
    CHECKMEM(d);

    if(wdim == 0) {
        // w is empty: there is only one IdSet, which contains
        // the samples, and goes in d[n-1] = d[nsamples].
        for(i=0; i < n-1; ++i)
            d[i] = IdSetSet_new(0);
        d[n-1] = IdSetSet_new(1);
        {
            IdSet *id = IdSet_new(nsamples, sample, 1.0);
            IdSet_sanityCheck(id, __FILE__, __LINE__);
            status = IdSetSet_add(d[n-1], id);
            if(status)
                ERR(status, "bad return from IdSetSet_add");
        }
        assert(1 == IdSetSet_size(d[n-1]));
        return d;
    }

    // w is not empty

    // If nsamples>0, the initial vector(s) are empty.
    for(i=0; i<nsamples; ++i)
        d[i] = IdSetSet_new(0);

    // Allocate non-empty vectors
    for(i=nsamples; i<n; ++i) {
        int j = i - nsamples;
        assert(w[j]);
        m = IdSetSet_size(w[j]);
        d[i] = IdSetSet_new(m);
    }

    // Fill non-empty vectors
    for(i=0; i < wdim; ++i) {
        IdSetSet_rewind(w[i]);
        for(IdSet *ids=IdSetSet_next(w[i]);
            ids;
            ids = IdSetSet_next(w[i])) {

            // Freed old ids and returns pointer to new one
            ids = IdSet_addSamples(ids, nsamples, sample);
            
            int nIds = IdSet_nIds(ids);
            status = IdSetSet_add(d[nIds], ids);
            if(status)
                ERR(status, "bad return from IdSetSet_add");
        }
        IdSetSet_empty_shallow(w[i]);
    }
    return d;
}

// Return a vector of IdSetSet objects. The i'th object contains sets
// of i descendants. The IdSet objects are generated by combining all
// compatible pairs from w0 and w1, and then adding the tipId_t values
// from array sample. On return, w0 and w1 are empty sets.
IdSetSet **get_descendants2(int dim0, IdSetSet **w0,
                            int dim1, IdSetSet **w1,
                            int nsamples, tipId_t *sample, int *newdim) {
    int i, j, n;

    // Is waiting room 0 empty?
    while(dim0 > 0 && IdSetSet_size(w0[dim0-1]) == 0)
        dim0 -= 1;

    // Is waiting room 1 empty?
    while(dim1 > 0 && IdSetSet_size(w1[dim1-1]) == 0)
        dim1 -= 1;

    // If either waiting room is empty, call get_descendants1.
    if(dim0 == 0)
        return get_descendants1(dim1, w1, nsamples, sample, newdim);

    if(dim1 == 0)
        return get_descendants1(dim0, w0, nsamples, sample, newdim);

    /*
     * Segment has two children. The returned value (dvec) is an array
     * of IdSetSet objects. The i'th object contains IdSet objects,
     * each of which of which contains i descendants. Each IdSet is
     * the union of (a) an entry from w0, (b) an entry from w1, and
     * (c) the samples (if any) in "sample". Some entries in w0 may be
     * incompatible with some in w1, because they represent mutually
     * exclusive outcomes of the same migration event. These mutually
     * exclusive pairs are not used.
     */
    if(nsamples || dim0 || dim1) {
        n = 1;
        n += nsamples;
        if(dim0)
            n += dim0-1;
        if(dim1)
            n += dim1-1;
    }else
        n = 0;

    if(n==0)
        return NULL;

    // Use PtrLst to begin with, because we can't predict how many
    // sets will be in each IdSetSet object. This avoids expensive
    // reallocations of the hash table.
    PtrLst *d[n];
    for(i=0; i<n; ++i)
        d[i] = PtrLst_new();

    IdSet *id0, *id1, *newid;

    for(i=0; i < dim0; ++i) {
        IdSetSet_rewind(w0[i]);
        for(id0=IdSetSet_next(w0[i]);
            id0;
            id0=IdSetSet_next(w0[i])) {
            
            IdSet_sanityCheck(id0, __FILE__, __LINE__);
            assert(IdSet_nIds(id0) == i);

            for(j=0; j < dim1; ++j) {
                IdSetSet_rewind(w1[j]);

                for(id1=IdSetSet_next(w1[j]);
                    id1;
                    id1=IdSetSet_next(w1[j])) {

                    IdSet_sanityCheck(id1, __FILE__, __LINE__);
                    assert(IdSet_nIds(id1) == j);

                    newid = IdSet_join(id0, id1, nsamples, sample);
                    if(newid == NULL)
                        continue;

                    IdSet_sanityCheck(newid, __FILE__, __LINE__);
                    
                    // non-NULL means id0 and id1 are compatible
                    // and can be added to the list of
                    // descendants.
                    int k = i+j+ nsamples;
                    assert(k == newid->nIds);
                    assert(k < n);
                    PtrLst_push(d[k], newid);
                }
            }
        }
        IdSetSet_empty_deep(w0[i]);
    }

    // Empty w1; w0 is already empty.
    for(j=0; j < dim1; ++j) {
        IdSetSet_empty_deep(w1[j]);
    }

    // Get rid of any empty lists at the top of d and reset n.
    while(n > 0 && 0 == PtrLst_length(d[n-1])) {
        PtrLst_free(d[n-1]);
        n -= 1;
    }
    *newdim = n;

    if(n == 0)
        return NULL;

    IdSetSet **dss = malloc(n * sizeof(dss[0]));

    // Convert lists of descendants into IdSetSet objects and
    // install in dss. Also free the entries of d[i].
    for(i=0; i<n; ++i) {
        int m = PtrLst_length(d[i]);
        dss[i] = IdSetSet_new(m);
        for(IdSet *id = PtrLst_pop(d[i]);
            id;
            id = PtrLst_pop(d[i])) {

            int status = IdSetSet_add(dss[i], id);
            if(status)
                ERR(status, "bad return from IdSetSet_add");
        }
        PtrLst_free(d[i]);
    }

    return dss;
}

static void coalescent_interval_length(int n, long double elen[n],
                                       long double eig[n-1],
                                       long double v) {
    // Calculate expected length, within the segment, of
    // coalescent intervals with 2,3,...,n lineages.  elen[i] is
    // expected length of interval with i+1 lineages.
    // elen[0] not set.
    MatCoal_ciLen(n-1, elen+1, eig);

    // Calculate expected length, elen[0], of interval with one
    // lineage.  elen[i] refers to interval with i+1 lineages.
    long double sum = 0.0;
    for(int i = n-1; i > 0; --i)
        sum += elen[i];
    elen[0] = v - sum;
}

static void project(int n, long double pr[n], long double eig[n-1]) {

    // Calculate prob of 2,3,...,n ancestors at ancent end of
    // segment. On return, pr[1] = prob[2], pr[n-1] = prob[n],
    // pr[0] not set.
    MatCoal_project(n-1, pr+1, eig);

    // Calculate probability of 1 line of descent. I'm doing
    // the sum in reverse order on the assumption that higher
    // indices will often have smaller probabilities.  pr[i]
    // is prob of i+1 lineages; pr[0] not set.
    long double sum=0.0;
    for(int i = n-1; i > 0; --i)
        sum += pr[i];
    pr[0] = 1.0L-sum;
}

/// Return 1 if each PtrLst in array w is empty; return 0 otherwise.
static int w_isempty(int dim, IdSetSet **w) {
    for(int i=0; i<dim; ++i) {
        if(IdSetSet_size(w[i]) > 0)
            return 0;
    }
    return 1;
}

/// Remove all references to samples from tree of populations.
/// Sets "visited" to 0 in every node.
void Segment_clear(Segment * self) {
    assert(self);
    for(int i = 0; i < self->nchildren; ++i) {
        assert(self->child[i]);
        Segment_clear(self->child[i]);
    }

    self->nsamples = 0;
    self->visited = 0;
    Segment_sanityCheck(self, __FILE__, __LINE__);
}

/// Return 1 if PopNode tree is empty of samples
int Segment_isClear(const Segment * self) {
    if(self == NULL)
        return 1;
    if(self->nsamples > 0)
        return 0;

    for(int i = 0; i < self->nchildren; ++i) {
        if(!Segment_isClear(self->child[i]))
            return 0;
    }
    return 1;
}

// Move all IdSet objects to parent ipar. Empties each list in array
// a.
static void mv_idsets_to_parent(Segment *self, int ipar, PtrLst **a) {
    assert(ipar < self->nparents);
    int iself = self_ndx(self, self->parent[ipar]);
    int status;

#ifndef NDEBUG    
    if(verbosity > 0) {
        PtrLst_rewind(a[1]);
        fprintf(stderr,"%s debug:", __func__);
        IdSet *set;
        while( (set = PtrLst_next(a[1])) != NULL) {
            IdSet_print(set, stderr);
        }
        putc('\n', stderr);
    }
#endif    

    assert(self == self->parent[ipar]->child[iself]);

    self->parent[ipar]->wdim[iself] = self->dim;
    for(int i=0; i < self->dim; ++i) {

        // parental waiting room
        IdSetSet *w = self->parent[ipar]->w[iself][i];

#ifndef NDEBUG
        IdSetSet_sanityCheck(w, __FILE__,__LINE__);
#endif
        
        int m = PtrLst_length(a[i]);
        IdSetSet_reserve(w, m);

        IdSet *id;
        while( (id = PtrLst_pop(a[i])) != NULL ) {
            status = IdSetSet_add(w, id);
            if(status)
                ERR(status, "bad return from IdSetSet_add");
        }
    }
}

/**
 * Move IdSet objects from list "src" into the relevant waiting room
 * of parent ipar. "nlin" is the number of lineages in each IdSet
 * object within "src". It also serves as the index into the parental
 * waiting room. On return, "src" is empty, and
 * self->parent[ipar]->w[iself][nlin] contains the IdSet objects that
 * were originally in "src". Here, "iself" is the index of "self" in
 * the parent's list of children.
 */
static void mv_to_waiting_room(Segment *self, PtrLst *src, int ipar,
                               int nlin) {
    int iself = self_ndx(self, self->parent[ipar]);
    assert(nlin < self->parent[ipar]->wdim[iself]);

    // parental waiting room
    IdSetSet *w = self->parent[ipar]->w[iself][nlin];

    int m = PtrLst_length(src);
    IdSetSet_reserve(w, m);

    IdSet *id;
    while( (id = PtrLst_pop(src)) != NULL ) {
        int status = IdSetSet_add(w, id);
        if(status)
            ERR(status, "bad return from IdSetSet_add");
    }
}

void Segment_print_d(Segment *self, const char *func, int line) {
    fprintf(stderr,"%s:%d: SEGNUM %d. d:\n",
            func,line, self->segnum);
    for(int i=0; i < self->dim; ++i) {
        fprintf(stderr, "%d@", i);
        IdSetSet_rewind(self->d[i]);
        for(IdSet *id = IdSetSet_next(self->d[i]);
            id;
            id = IdSetSet_next(self->d[i])) {

            IdSet_print(id, stderr);

        }
        putc('\n', stderr);
    }
}

/// This sorts tipId_t values into numerical order, unlike the
/// more complex comparison function, compare_tipId, which is defined
/// in lblndx.c.
static int tipidcmp(const void *vx, const void *vy) {
    tipId_t const * x = vx;
    tipId_t const * y = vy;
    if(*x > *y)
        return 1;
    if(*x < *y)
        return -1;
    return 0;
}
