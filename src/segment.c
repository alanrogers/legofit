/**
 * @file segment.c
 * @author Alan R. Rogers
 * @brief A single segment of a population tree.
 *
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "binary.h"
#include "branchtab.h"
#include "comb.h"
#include "error.h"
#include "idset.h"
#include "matcoal.h"
#include "misc.h"
#include "parstore.h"
#include "ptrlst.h"
#include "ptrptrmap.h"
#include "ptrvec.h"
#include "segment.h"
#include "setpart.h"
#include <stdlib.h>
#include <string.h>

// The total number of samples.
//static int total_samples = 0;

typedef struct MigDat MigDat;
typedef struct CombDat CombDat;
typedef struct SetPartDat SetPartDat;

// For assigning initial tipId_t values within segments.
//static tipId_t currTipId = 0;

// Data manipulated by visitComb function
struct CombDat {
    double contrib; // Pr[site pattern]*E[len of interval]
    PtrVec *d;
    BranchTab *branchtab;
    int dosing;    // do singleton site patterns if nonzero
};

// Data manipulated by migrate function
struct MigDat {
    int nMigrants, nNatives;
    double pr; // probability of 
    PtrLst *migrants, *natives;
    PtrVec *a;
    unsigned migrationEvent, outcome;
};

// Data manipulated by visitSetPart function.
struct SetPartDat {
    unsigned nparts; // Number of parts in partition

    double prior; // prob of k ancestors given n descendants
    long double lnconst;  // log of constant in Durrett's theorem 1.5
    double elen;     // E[len of interval]

    // a is a list of IdSet objects for sets of ancestors
    PtrLst *a;

    // PtrVec_get(d, i) is the IdSet for the i'th set of descendants
    PtrVec *d;

    BranchTab *branchtab;
    int dosing;      // do singleton site patterns if nonzero
};

// One segment of a population network. This version works
// with MCTree.
struct Segment {
    int            nparents, nchildren, nsamples;
    int            visited;     // for traversal algorithm
    double         twoN;        // ptr to current pop size
    double         start, end;  // duration of this Segment
    double         mix;         // ptr to frac of pop derived from parent[1]

    // indices into ParStore array
    int twoN_i, start_i, end_i, mix_i;

    struct Segment *parent[2];
    struct Segment *child[2];

    // max number of lineages in this segment.
    int max;

    tipId_t sample[MAXSAMP]; // array of length nsamples

    // Waiting rooms.  Each child loads IdSet objects into one of two
    // waiting rooms. w[i] is the i'th waiting room, where i is in
    // [0,1]. w[i] is the waiting room for child i. The number of
    // waiting rooms equals the number of children. w[i][j] is a list
    // of IdSet objects, each of which contains i+1 tipId_t values.
    // wmax[i] is the maximum number of lineages in the i'th waiting
    // room, so 0 <= j < wmax[i]. If there are no waiting rooms
    // (i.e. if nchildren=0), then the descendants at the beginning of the
    // segment consist only of those in array "sample". If there is
    // one waiting room, then the ids in "sample" are added to each of
    // the sets in w[0], and the resulting list of IdSet objects
    // represents the state at the beginning of the segment. If there
    // are two waiting rooms, then we begin with all pairs of IdSet
    // objects in w[0] and w[1]. To each pair, we add the ids in
    // "samples" to obtain the list of IdSet objects at the beginning
    // of the segment. 
    int wmax[2];
    PtrLst *w[2][MAXSAMP];

    // d[i] is a pointer to a PtrVec, which is allocated within
    // Segment_coalesce. This PtrVec holds IdSet objects that contain
    // i+1 tipId_t values, each representing a descendant at the
    // recent end of the Segment. The number of PtrVec
    // pointers--i.e. the dimension of array d--is self->max.
    PtrVec **d;
};

static PtrVec **get_descendants1(int dim, PtrLst **w, int nsamples,
                                 tipId_t *sample, int *newdim);
static PtrVec **get_descendants2(int dim0, PtrLst **w0,
                                 int dim1, PtrLst **w1,
                                 int nsamples, tipId_t *sample, int *newdim);
static void mv_idsets_to_parent(Segment *self, int ipar, PtrLst **a);
static void coalescent_interval_length(int n, double elen[n],
                                       double eig[n], double v);
static void project(int n, double pr[n], double eig[n]);
int    visitComb(int d, int ndx[d], void *data);
int    visitSetPart(unsigned n, unsigned a[n], void *data);
int    visitMig(int nmig, int *migndx, void *data);
static void  unlink_child(Segment *child, Segment *parent);
static int Segment_coalesceFinite(Segment *self, double v, int dosing,
                                  BranchTab *branchtab);
static int Segment_coalesceInfinite(Segment *self, double v, int dosing,
                                    BranchTab *branchtab);
static void Segment_duplicate_nodes(Segment *old, PtrPtrMap *ppm);
static int Segment_equals_r(Segment *a, Segment *b);
static int self_ndx(Segment *self, Segment *parent);
PtrVec *PtrVec_from_PtrLst(PtrVec *to, PtrLst *from);
static int w_isempty(int dim, PtrLst **w);

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

    self->twoN_i = twoN_i;
    self->start_i = start_i;
    self->end_i = -1;
    self->mix_i = -1;

    self->twoN = ParStore_getVal(ps, twoN_i);
    self->start = ParStore_getVal(ps, start_i);
    self->end = INFINITY;
    self->mix = 0.0;

    for(int i=0; i < MAXSAMP; ++i) {
        self->w[0][i] = PtrLst_new();
        self->w[1][i] = PtrLst_new();
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

    for(int i=0; i < MAXSAMP; ++i) {
        PtrLst_free(self->w[0][i]);
        PtrLst_free(self->w[1][i]);
    }

    free(self);
}

/// Traverse network, removing segments with no children and
/// no samples. These do not contribute to the calculation.
void Segment_prune(Segment *self) {
    if(self->nchildren > 1)
        Segment_prune(self->child[1]);
    if(self->nchildren > 0)
        Segment_prune(self->child[0]);

    if(self->nchildren == 0 && self->nsamples == 0)
        Segment_free(self);
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
    // and identify the root of the duplicated network.
    for(unsigned i=0; i < nnodes; ++i) {
        Segment *old = old_nodes[i];
        Segment *new = PtrPtrMap_get(ppm, old, &status);
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
    CombDat *dat = (CombDat *) data;


    unsigned nIdSets = PtrVec_length(dat->d);

    assert(nIdSets > 0);
    
    for(unsigned i=0; i < nIdSets; ++i) {
        IdSet *ids = PtrVec_get(dat->d, i);
        
        // sitepat is the union of the current set of descendants, as
        // described in ndx.
        tipId_t sitepat = 0;
        for(int j=0; j < d; ++j) {
            assert(ndx[j] < ids->nIds);
            assert(ndx[j] >= 0);
#ifndef NDEBUG
            if(ids->tid[ndx[j]] == 0) {
                fprintf(stderr,"%s:%d: j=%d d=%d nIds=%d ids->tid[%u] = 0\n",
                        __FILE__,__LINE__, j, d, ids->nIds, ndx[j]);
            }
#endif            
            assert(ids->tid[ndx[j]]);
            sitepat |= ids->tid[ndx[j]];
        }
        
        assert(sitepat > 0);

        // Skip singletons unless data->dosing is nonzero
        if(!dat->dosing && isPow2(sitepat))
            continue;

      
        // Increment BranchTab entry for current sitepat value.
        BranchTab_add(dat->branchtab, sitepat,
                      ids->p * dat->contrib);
    }
    return 0;
}

/// Visit a set partition. n is the number of descendants, a[i] is the
/// index of the ancestor of the i'th descendant.
int visitSetPart(unsigned n, unsigned a[n], void *data) {
    SetPartDat *vdat = (SetPartDat *) data;
    int status=0;

    // Calculate the size, c[i], of each part. In other words,
    // the number of descendants of each ancestor.
    unsigned k=vdat->nparts, c[k];
    memset(c, 0, k*sizeof(c[0]));
    for(int i=0; i<n; ++i) {
        assert(a[i] < k);
        ++c[a[i]];
    }

    double p = probPartition(k, c, vdat->lnconst);
    p *= vdat->prior;

    // nIds is the number of IdSet objects, each representing
    // a set of n descendants.
    unsigned nIds = PtrVec_length(vdat->d);
    for(unsigned i=0; i < nIds; ++i) {
        IdSet *descendants = PtrVec_get(vdat->d, i);
        tipId_t sitepat[k];
        memset(sitepat, 0, k*sizeof(tipId_t));

        assert(n == descendants->nIds);

        // Loop over descendants, creating a sitepat for each
        // ancestor. a[i] is the index of the ancestor of the i'th
        // descendant. sitepat[j] is the site pattern of the j'th
        // ancestor.
        for(int j=0; j<n; ++j)
            sitepat[a[j]] |= descendants->tid[j];

        // Loop over ancestors, i.e. over site patterns, adding
        // to the corresponding entry in BranchTab.
        for(int j=0; j<k; ++j)
            BranchTab_add(vdat->branchtab, sitepat[j],
                          p * vdat->elen * descendants->p);

        // Add the current set partition to the list of ancestral
        // states.

        IdSet *ancestors = IdSet_new(k, sitepat, p * descendants->p);
        IdSet_sanityCheck(ancestors, __FILE__, __LINE__);
        IdSet_copyMigOutcome(ancestors, descendants);
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
/// indices of the current set of migrants. Its length is nmig.
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

    // Arrays of tipId_values also have an extra entry
    // to guard against problems of zero length.
    tipId_t migid[1 + nmig], natid[1 + nnat];

    // number of sets of ancestors
    int nSets = PtrVec_length(mdat->a);

    for(int i_set=0; i_set < nSets; ++i_set) {
        IdSet *set = PtrVec_get(mdat->a, i_set);

        assert(IdSet_nIds(set) == nmig + nnat);

        for(i=0; i<nmig; ++i)
            migid[i] = set->tid[migndx[i]];

        if(nmig == 0) {
            nmig = 1;
            migid[0] = 0; // empty set
        }

        for(i=0; i < nnat; ++i)
            natid[i] = set->tid[natndx[i]];

        if(nnat == 0) {
            nnat = 1;
            natid[0] = 0; // empty set
        }

        // Create IdSet objects for migrants and natives
        IdSet *mig = IdSet_new(nmig, migid, set->p);
        IdSet_addMigEvent(mig, mdat->migrationEvent,
                          mdat->outcome, mdat->pr);
        IdSet_sanityCheck(mig, __FILE__, __LINE__);
        PtrLst_push(mdat->migrants, mig);

        IdSet *nat = IdSet_new(nnat, natid, set->p);
        IdSet_addMigEvent(nat, mdat->migrationEvent,
                          mdat->outcome, mdat->pr);
        IdSet_sanityCheck(nat, __FILE__, __LINE__);
        PtrLst_push(mdat->natives, nat);
    }
    mdat->outcome += 1;

    return 0;
}


/// If to==NULL, then return a new PtrVec, containing all the
/// pointers in "from". On return, "from" is empty, but the
/// PtrLst object itself is not freed.
///
/// If to!=NULL, then resize "to" so that it is large enough to hold
/// all the pointers in "from", and move all pointers from "from" into
/// "to". Return "to".
PtrVec *PtrVec_from_PtrLst(PtrVec *to, PtrLst *from) {
    unsigned nIdSets = PtrLst_length(from);

    if(to) {
        // resize "to"
        assert(PtrVec_length(to) == 0);
        if(PtrVec_resize(to, nIdSets)) {
            fprintf(stderr,"%s:%d: allocation error\n",
                    __FILE__,__LINE__);
            exit(EXIT_FAILURE);
        }
    }else {
        // fresh allocation
        to = PtrVec_new(nIdSets);
    }
    
    for(IdSet *idset = PtrLst_pop(from); idset; idset = PtrLst_pop(from))
        PtrVec_push(to, idset);

    return to;
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

    for(int i=0; i < MAXSAMP; ++i) {
        new->w[0][i] = PtrLst_new();
        new->w[1][i] = PtrLst_new();
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

static int Segment_coalesceFinite(Segment *self, double v, int dosing,
                                  BranchTab *branchtab) {
    assert(self->max > 0);
    double pr[self->max], elen[self->max];
    int n, i, k, iself, status=0;

    // Array of lists of ancestors. a[k-1] is the list for sets of k
    // ancestors.  
    PtrLst *a[self->max];
    for(i=0; i < self->max; ++i)
        a[i] = PtrLst_new();
    
    SetPartDat sd = {.branchtab = branchtab,
                     .dosing = dosing,
    };

    // Handle case of a segment with a single lineage
    if(PtrVec_length(self->d[0]) > 0) {

        int nIds = PtrVec_length(self->d[0]);
        for(i=0; i < nIds; ++i) {
            IdSet *d = PtrVec_get(self->d[0], i);
            assert(d->nIds == 1);
            BranchTab_add(branchtab, d->tid[0], d->p * v);

            IdSet *new = IdSet_dup(d);
            PtrLst_push(a[0], new);
        }
    }
    
    // Outer loop over numbers of descendants.
    // Calculate probabilities and expected values, p[1] and elen.
    for(n=2; n <= self->max; ++n) {

        // eigenvalues of transient states
        double eig[n-1];
        MatCoal_eigenvals(n-1, eig, v);

        // Calculate the expected length, elen[i], of the subinterval
        // containing i+1 lineages.
        coalescent_interval_length(n, elen, eig, v);
        
        // Calculate pr[i], the probability of i+1 lineages at the
        // ancient end of the segment.
        project(n, pr, eig);

        sd.d = self->d[n-1];

        // Loop over number, k, of ancestors.
        // Include k=1, because this is a finite segment.
        for(k=1; k <= n; ++k) {
            sd.a = a[k-1];
            assert(0 == PtrLst_length(sd.a));
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
        mv_idsets_to_parent(self, 0, a);
    }else{
        assert(self->nparents == 2);
        assert(self->mix_i >= 0);
        assert(self->mix >= 0.0);
        assert(self->mix <= 1.0);

        // natives go to parent[0], migrants to parent[1]

        // Set dimension of w arrays in the two parents
        for(int ipar = 0; ipar < 2; ++ipar) {
            iself = self_ndx(self, self->parent[ipar]);
            self->parent[ipar]->wmax[iself] = self->max;
        }

        //  P[x] = (k choose x) * m^x * (1-m)^(k-x)
        MigDat msd = {
                      .migrants = PtrLst_new(),
                      .natives = PtrLst_new(),
                      .a = NULL,
                      .migrationEvent = nextMigrationEvent(),
                      .outcome = 0
        };
        for(k=1; k <= self->max; ++k) {
            msd.a = PtrVec_from_PtrLst(msd.a, a[k-1]); // empties a[k-1]
            for(long x=0; x <= k; ++x) {
                // prob that x of k lineages are migrants
                long double lnpr = lbinom(k, x)
                    + x*logl(self->mix)
                    + (k-x)*logl(1.0-self->mix);
                
                msd.pr = expl(lnpr);
                msd.nMigrants = x;
                msd.nNatives = k - x;
                status = traverseComb(k, x, visitMig, &msd);
                if(status)
                    return status;

                // transfer natives
                if(msd.nNatives > 0) {
                    iself = self_ndx(self, self->parent[0]);
                    assert(k-x-1 < self->parent[0]->wmax[iself]);
                    PtrLst_move(self->parent[0]->w[iself][k-x-1], msd.natives);
                }

                // transfer migrants
                if(msd.nMigrants > 0) {
                    iself = self_ndx(self, self->parent[1]);
                    assert(x > 0);
                    assert(x-1 < self->parent[1]->wmax[iself]);
                    PtrLst_move(self->parent[1]->w[iself][x-1], msd.migrants);
                }
            }
        }
        PtrLst_free(msd.migrants);
        PtrLst_free(msd.natives);
        for(IdSet *s=PtrVec_pop(msd.a); s; s=PtrVec_pop(msd.a))
            IdSet_free(s);
        PtrVec_free(msd.a);
    }

    for(i=0; i < self->max; ++i)
        PtrLst_free(a[i]);

    return status;
}

static int Segment_coalesceInfinite(Segment *self, double v, int dosing,
                                    BranchTab *branchtab) {
    assert(self->max > 0);
    double elen[self->max];
    int n, status=0;

    CombDat cd = {.branchtab = branchtab,
                  .dosing = dosing
    };

    // Outer loop over numbers of descendants.
    // Calculate expected branch lengths, elen.
    for(n=2; n <= self->max; ++n) {

        // In an infinite interval, eigenvalues of all transient
        // states are zero.
        double eig[n-1];
        memset(eig, 0, (n-1)*sizeof(eig[0]));
        
        // Calculate the expected length, elen[i], of the subinterval
        // containing i+1 lineages.
        coalescent_interval_length(n, elen, eig, v);

        cd.d = self->d[n-1];

        // Loop over number, k, of ancestors.
        // Exclude k=1, because this is an infinite segment.
        for(int k=2; k <= n; ++k) {
            
            // portion of log Qdk that doesn't involve d
            long double lnconst = logl(k) - lbinom(n-1, k-1);

            // Within each interval, there can be ancestors
            // with 1 descendant, 2, 3, ..., n-k+2.
            for(int d=1; d <= n-k+1; ++d) {
                long double lnprob = lnconst
                    + lbinom(n-d-1, k-1) - lbinom(n,d);

                // probability of site pattern
                cd.contrib = (double) expl(lnprob);

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
    int status=0;

    if(self->visited)
        return 0;
    self->visited = 1;

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

    // Set self->d, the array of descendants, and self->max, the
    // dimension of this array.
    switch(self->nchildren) {
    case 0:
        self->d = get_descendants1(0, NULL,
                                   self->nsamples, self->sample,
                                   &self->max);
        assert(self->max > 0);
        break;
    case 1:
        self->d = get_descendants1(self->wmax[0], self->w[0],
                                   self->nsamples, self->sample,
                                   &self->max);
        assert(self->max > 0);
        assert(w_isempty(self->wmax[0], self->w[0]));
        break;
    case 2:
        self->d = get_descendants2(self->wmax[0], self->w[0],
                                   self->wmax[1], self->w[1],
                                   self->nsamples, self->sample,
                                   &self->max);
        assert(self->max > 0);
        assert(w_isempty(self->wmax[0], self->w[0]));
        assert(w_isempty(self->wmax[1], self->w[1]));
        break;
    default:
        fprintf(stderr,"%s:%d: illegal number of children: %d\n",
                __FILE__,__LINE__, self->nchildren);
        exit(EXIT_FAILURE);
    }
    assert(self->max > 0);

    fprintf(stderr,"%s: wmax=[%d,%d] max=%d 2N=%lg t=[%lg,%lg]\n",
            __func__, self->wmax[0], self->wmax[1], self->max,
            self->twoN, self->start, self->end);

    for(int i=0; i < self->max; ++i) {
        fprintf(stderr,"%s:%d: nIds=%d\n", __FILE__,__LINE__, i+1);
        for(int j=0; j < PtrVec_length(self->d[i]); ++j)
            IdSet_print(PtrVec_get(self->d[i], j), stderr);
    }

    if(self->end_i == -1) {
        status = Segment_coalesceInfinite(self, INFINITY, dosing, branchtab);
    }else{
        double v = self->end - self->start;
        status = Segment_coalesceFinite(self, v, dosing, branchtab);
    }

    // Free IdSet objects of descendants
    for(int i=0; i < self->max; ++i) {
        for(IdSet *ids=PtrVec_pop(self->d[i]);
            ids;
            ids = PtrVec_pop(self->d[i])) {
            
            IdSet_free(ids);

        }
        PtrVec_free(self->d[i]);
    }
    free(self->d);
    self->d = NULL;

    return status;
}

/// Return a newly-allocated array of PtrVec vectors. The i'th
/// vector contains all the IdSet entries in w[i], with each entry
/// augmented by the entries (if any) in vector "sample". If there are
/// no entries in any of the w[i], the returned array contains a
/// single non-empty entry, with the entries of "samples". If "w" and
/// "samples" are both empty, *newdim will equal 0 and the function
/// returns NULL. Otherwise, *newdim is the dimension of the returned
/// array.
///
/// On entry, dim is the dimension of array w, and nsamples is the
/// dimension of array sample. On return *newdim is the dimension of
/// the newly-allocated array returned by the function.
static PtrVec **get_descendants1(int dim, PtrLst **w, int nsamples,
                                 tipId_t *sample, int *newdim) {
    int i, n, m;

    while(dim > 0 && PtrLst_length(w[dim-1]) == 0)
        dim -= 1;

    n = dim + nsamples;
    *newdim = n;
    if(n==0) {
        return NULL;
    }

    PtrVec **d = malloc(n * sizeof(PtrVec *));
    CHECKMEM(d);

    if(dim == 0) {
        // w is empty
        for(i=0; i < n-1; ++i)
            d[i] = PtrVec_new(0);
        d[n-1] = PtrVec_new(1);
        {
            IdSet *id = IdSet_new(n, sample, 1.0);
            IdSet_sanityCheck(id, __FILE__, __LINE__);
            PtrVec_push(d[n-1], id);
        }
        assert(1 == PtrVec_length(d[n-1]));
        return d;
    }

    // w is not empty

    // If nsamples>0, the initial vector(s) are empty.
    for(i=0; i<nsamples; ++i)
        d[i] = PtrVec_new(0);

    // Allocate non-empty vectors
    for(i=nsamples; i<n; ++i) {
        int j = i - nsamples;
        assert(w[j]);
        m = PtrLst_length(w[j]);
        d[i] = PtrVec_new(m);
    }

    // Fill non-empty vectors
    for(i=0; i < dim; ++i) {
        IdSet *ids = PtrLst_pop(w[i]);
        while(ids) {
            ids = IdSet_addSamples(ids, nsamples, sample);
            int nIds = IdSet_nIds(ids);
            PtrVec_push(d[nIds-1], ids);
            ids = PtrLst_pop(w[i]);
        }
    }
    return d;
}

// Return a vector of vectors of IdSet objects. The i'th vector
// contains sets of i+1 descendants. The IdSet objects are generated
// by combining all compatible pairs from w0 and w1, and then adding
// the tipId_t values from array sample. On return, w0 and w1 are
// arrays of empty lists.
PtrVec **get_descendants2(int dim0, PtrLst **w0,
                          int dim1, PtrLst **w1,
                          int nsamples, tipId_t *sample, int *newdim) {
    int i, j, n;

    // Is waiting room 0 empty?
    while(dim0 > 0 && PtrLst_length(w0[dim0-1]) == 0)
        dim0 -= 1;

    // Is waiting room 1 empty?
    while(dim1 > 0 && PtrLst_length(w1[dim1-1]) == 0)
        dim1 -= 1;

    // If either waiting room is empty, call get_descendants1.
    if(dim0 == 0)
        return get_descendants1(dim1, w1, nsamples, sample, newdim);

    if(dim1 == 0)
        return get_descendants1(dim0, w0, nsamples, sample, newdim);

    /*
     * Segment has two children. The returned value (dvec) is an array
     * of vectors. The i'th vector contains sets of i+1
     * descendants. Each such set is the union of (a) an entry from
     * w0, (b) an entry from w1, and (c) the samples (if any) in
     * "sample". Some entries in w0 may be incompatible with some in
     * w1, because they represent mutually exclusive outcomes of the
     * same migration event. These mutually exclusive pairs are not
     * used. For this reason, we cannot figure out in advance how many
     * IdSet objects will be in each entry of dvec[i]. To solve this
     * problem, the algorithm first fills an array of linked
     * lists. Then each linked list is converted into an array to form
     * an entry of dvec.
     */
    n = dim0 + dim1 + nsamples;
    assert(n > 0);

    PtrLst *d[n];
    for(i=0; i<n; ++i)
        d[i] = PtrLst_new();

    IdSet *id0, *id1, *newid;

    for(i=0; i < dim0; ++i) {
        for(id0=PtrLst_pop(w0[i]);
            id0;
            id0=PtrLst_pop(w0[i])) {
            
            IdSet_sanityCheck(id0, __FILE__, __LINE__);

            for(j=0; j < dim1; ++j) {
                PtrLst_rewind(w1[j]);

                for(id1=PtrLst_next(w1[j]);
                    id1;
                    id1=PtrLst_next(w1[j])) {

                    IdSet_sanityCheck(id1, __FILE__, __LINE__);
                    newid = IdSet_join(id0, id1, nsamples, sample);
                    if(newid == NULL)
                        continue;
                    IdSet_sanityCheck(newid, __FILE__, __LINE__);
                    
                    // non-NULL means id0 and id1 are compatible
                    // and can be added to the list of
                    // descendants.
                    int k = i+j+ nsamples + 2;
                    assert(k == newid->nIds);
                    PtrLst_push(d[k-1], newid);
                }
            }
            IdSet_free(id0);
        }
    }

    // Get rid of any empty lists at the top of d and reset n.
    while(n > 0 && 0 == PtrLst_length(d[n-1])) {
        PtrLst_free(d[n-1]);
        n -= 1;
    }
    *newdim = n;

    if(n == 0)
        return NULL;

    PtrVec **dvec = malloc(n * sizeof(PtrVec *));

    // Convert lists of descendants into vectors and
    // install in dvec. Also free the entries of d[i].
    for(i=0; i<n; ++i) {
        dvec[i] = PtrVec_from_PtrLst(NULL, d[i]);
        PtrLst_free(d[i]);
    }

    // Empty w1. w0 is already empty.
    for(i=0; i < dim1; ++i) {
        for(id1=PtrLst_pop(w1[i]);
            id1;
            id1=PtrLst_pop(w1[i])) {

            IdSet_free(id1);
        }
    }

    return dvec;
}

static void coalescent_interval_length(int n, double elen[n],
                                       double eig[n], double v) {
    // Calculate expected length, within the segment, of
    // coalescent intervals with 2,3,...,n lineages.  elen[i] is
    // expected length of interval with i-1 lineages.
    // elen[0] not set.
    MatCoal_ciLen(n-1, elen+1, eig);

    // Calculate expected length, elen[0], of interval with one
    // lineage.  elen[i] refers to interval with i+1 lineages.
    double sum = 0.0;
    for(int i=n; i >= 2; --i)
        sum += elen[i-1];
    elen[0] = v - sum;

}

static void project(int n, double pr[n], double eig[n]) {
    // Calculate prob of 2,3,...,n ancestors at ancent end of
    // segment. On return, pr[1] = prob[2], pr[n-1] = prob[n],
    // pr[0] not set.
    MatCoal_project(n-1, pr+1, eig);

    // Calculate probability of 1 line of descent. I'm doing
    // the sum in reverse order on the assumption that higher
    // indices will often have smaller probabilities.  pr[i]
    // is prob of i+1 lineages; pr[0] not set.
    double sum=0.0;
    for(int i=n; i > 1; --i)
        sum += pr[i-1];
    pr[0] = 1.0-sum;

}

/// Return 1 if each PtrLst in array w is empty; return 0 otherwise.
static int w_isempty(int dim, PtrLst **w) {
    for(int i=0; i<dim; ++i) {
        if(PtrLst_length(w[i]) > 0)
            return 0;
    }
    return 1;
}

// Move all IdSet objects to parent ipar. Empties each list
// in array a.
static void mv_idsets_to_parent(Segment *self, int ipar, PtrLst **a) {
    int iself = self_ndx(self, self->parent[ipar]);

    assert(self == self->parent[ipar]->child[iself]);

    self->parent[ipar]->wmax[iself] = self->max;
    for(int i=0; i < self->max; ++i) {
        PtrLst_move(self->parent[ipar]->w[iself][i], a[i]);
    }
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
    //memset(self->sample, 0, sizeof(self->sample));
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
