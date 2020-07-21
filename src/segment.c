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
#include "setpart.h"
#include "matcoal.h"
#include "misc.h"
#include "segment.h"
#include <stdlib.h>
#include <string.h>

// The total number of samples.
static int total_samples = 0;

typedef struct MigDat MigDat;
typedef struct CombDat CombDat;
typedef struct SetPartDat SetPartDat;

// For assigning initial tipId_t values within segments.
static tipId_t currTipId = 0;

// Data manipulated by visitComb function
struct CombDat {
    double contrib; // Pr[site pattern]*E[len of interval]
    PtrVec *d;
    BranchTab *branchtab;
    int dosing;    // do singleton site patterns if nonzero
};

// Data manipulated by migrate function
struct MigDat {
    double pr; // probability
    PtrLst *migrants, *natives;
    PtrVec *a;
    unsigned migrationEvent;
};

// Data manipulated by visitSetPart function.
struct SetPartDat {
    unsigned nparts; // Number of parts in partition
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
    double         twoN;        // ptr to current pop size
    double         start, end;  // duration of this PopNode
    double         mix;         // ptr to frac of pop derived from parent[1]

    // indices into ParStore array
    int twoN_i, start_i, end_i, mix_i;

    struct Segment *parent[2];
    struct Segment *child[2];

    tipId_t sample[MAXSAMP]; // array of length nsamples

    int max;       // max number of lineages in segment

    // d[i] is a vector sets of i+1 descendants a[i] is a vector of
    // sets of i+1 ancestors Array d has dimension "max", which is not
    // known until the segments are linked into a
    // network. Consequently, the constructor sets its values to NULL,
    // and it is allocated only after the network has been assembled.
    PtrVec         **d;

    // Waiting rooms.  Each child loads IdSet objects into one of two
    // waiting rooms.  The descendants at the beginning of the segment
    // then consist of all pairs formed by two IdSets, one from each
    // waiting room. nw is the number of waiting room. The first child
    // fills w[0] and increments nw.
    int nw;
    PtrVec **w[2];

    // p[0][i] is prob there are i+1 lineages at recent end of segment
    // p[1][i] is analogous prob for ancient end of interval.
    double p[2][MAXSAMP];
};

int    visitComb(int d, int ndx[d], void *data);
int    visitSetPart(unsigned n, unsigned a[n], void *data);
int    visitMig(unsigned n, unsigned a[n], void *data);
static int Segment_coalesceFinite(Segment *self, double v, int dosing,
                                  BranchTab *branchtab);
static int Segment_coalesceInfinite(Segment *self, double v, int dosing,
                                    BranchTab *branchtab);
PtrVec *PtrLst_to_PtrVec(PtrLst *old);

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

    // other fields are not initialized by Segment_new.

    return self;
}

double Segment_twoN(Segment *self) {
    return self->twoN;
}

double Segment_start(Segment *self) {
    return self->start;
}

double Segment_end(Segment *self) {
    return self->end;
}

double Segment_mix(Segment *self) {
    return self->mix;
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
    REQUIRE(self->max >= 0 && self->max < MAXSAMP, file, lineno);
    for(int i=0; i < self->max; ++i) {
        REQUIRE(self->p[0][i] >= 0.0, file, lineno);
        REQUIRE(self->p[1][i] >= 0.0, file, lineno);
        REQUIRE(self->p[0][i] <= 1.0, file, lineno);
        REQUIRE(self->p[1][i] <= 1.0, file, lineno);
        IdSet_sanityCheck(self->ids[0][i], file, lineno);
        IdSet_sanityCheck(self->ids[1][i], file, lineno);
    }
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

void     Segment_print(FILE * fp, void * vself, int indent) {
    Segment *self = vself;
    for(int i = 0; i < indent; ++i)
        fputs("   ", fp);
    fprintf(fp, "%p twoN=%lf ntrval=(%lf,", self, self->twoN, self->start);
    fprintf(fp, "%lf)\n", self->end);

    for(int i = 0; i < self->nchildren; ++i)
        Segment_print(fp, self->child[i], indent + 1);
}

// Call from root Segment to set values of "max" within each
// Segment in network, and to allocate arrays "a" and "d"
void Segment_allocArrays(Segment *self) {
    int status;

    // If this segment has already been set, return immediately.
    if(self->max > 0)
        return;

    // Set children before self
    if(self->nchildren > 0)
        Segment_allocArrays(self->child[0]);
    if(self->nchildren == 2)
        Segment_allocArrays(self->child[1]);

    /* Absent migration, the maximum number of samples in the current
     * segment would be the sum of nsamples plus the maxima of the
     * children. But with migration, some of the potential lineages in
     * child 0 may be the same as some of those in child 1,
     * representing different and mutually exclusive outcomes of the
     * migration process. Therefore the sum of the children's maxima
     * will sometimes overestimate the max of the current node. To
     * guard against this, we make sure that self->max does not exceed
     * the total number of samples in the entire data set. */
    self->max = self->nsamples;
    if(self->nchildren > 0)
        self->max += self->child[0]->max;
    if(self->nchildren == 2)
        self->max += self->child[1]->max;

    if(self->max > total_samples)
        self->max = total_samples;

    // This could happen if a user creates a tip segment with no
    // samples.
    if(self->max == 0)
        return;

    // max number of ancestors equals max of descendants, because
    // the ancestors = descendants if no coalescent events happen.
    self->a = malloc(self->max * sizeof(PtrVec *));
    CHECKMEM(self->a);
    self->d = malloc(self->max * sizeof(PtrVec *));
    CHECKMEM(self->d);

    // The number of sets of k (=i+1) ancestors is the number
    // of ways of partitioning descendants into k subsets, and
    // this is given by Stirling's number of the 2nd kind.
    for(int i=0; i < self->max; ++i) {
        // number of ways to partition max things into i+1 subsets
        long unsigned n = stirling2(max, i+1);
        self->d[i] = PtrVec_new(n);
        self->a[i] = PtrVec_new(n);
    }

    // Create IdSet objects for tip Segments and store them
    // in self->a.
    if(self->nchildren==0 && self->nsamples > 0) {
        IdSet *idset = IdSet_new(self->nsamples, self->samples, 1.0);
        status = PtrVec_push(self->a[self->nsamples-1], idset);
        if(status) {
            fprintf(stderr,"%s:%d: %s\n", __FILE__,__LINE__,
                    strerror(status));
            exit(EXIT_FAILURE);
        }
    }
}

static int Segment_coalesceFinite(Segment *self, double v, int dosing,
                                  BranchTab *branchtab) {
    double pr[self->max], elen[self->max];
    int n, i, status=0;
    IdSet *ids;

    // Array of lists of ancestors. a[k-1] is the list for sets of k
    // ancestors.  
    PtrLst *a[self->max];
    for(i=0; i < self->max; ++i)
        a[i] = PtrLst();
    
    SetPartDat sd = {.branchtab = branchtab,
                     .dosing = dosing,
    };

#ifndef NDEBUG
    for(i=0; i < self->max; ++i)
        assert(self->ids[1][i] == NULL);
#endif    

    memset(self->p[1], 0, MAXSAMP*sizeof(double));

    // If there is only one line of descent, no coalescent events
    // are possible, so p[1][0] is at least as large as p[0][0].
    self->p[1][0] = self->p[0][0];

    // Outer loop over numbers of descendants.
    // Calculate probabilities and expected values, p[1] and elen.
    for(n=2; n <= self->max; ++n) {

        // eigenvalues of transient states
        double eig[n-1];
        MatCoal_eigenvals(n-1, eig, v);
        
        // Calculate expected length, within the segment, of
        // coalescent intervals with 2,3,...,n lineages.  elen[i] is
        // expected length of interval with i-1 lineages.
        // elen[0] not set.
        MatCoal_ciLen(n-1, elen+1, eig);

        // Calculate expected length, elen[0], of interval with one
        // lineage.  elen[i] refers to interval with i+1 lineages.
        double sum = 0.0;
        for(i=n; i >= 2; --i)
            sum += elen[i-1];
        elen[0] = v - sum; // elen[0] is infinite if v is.

        // Multiply elen by probability that segment had n lineages on
        // recent end of segment.
        for(i=0; i < n; ++i)
            elen[i] *= self->p[0][n-1];

        // Calculate prob of 2,3,...,n ancestors at ancent end of
        // segment. On return, pr[1] = prob[2], pr[n-1] = prob[n],
        // pr[0] not set.
        MatCoal_project(n-1, pr+1, eig);

        // Calculate probability of 1 line of descent. I'm doing
        // the sum in reverse order on the assumption that higher
        // indices will often have smaller probabilities.  pr[i]
        // is prob of i+1 lineages; pr[0] not set
        sum=0.0;
        for(i=n; i > 1; --i)
            sum += pr[i-1];
        self->p[1][0] += self->p[0][n-1] * (1.0-sum);

        // Add probs of 2..n lines of descent
        // pr[i] is prob of i+2 lineages
        // self->p[1][i] is prob if i+1 lineages
        for(i=2; i <= n; ++i)
            self->p[1][i-1] += self->p[0][n-1] * pr[i-1];

        sd.d = self->d[n-1];

        // Loop over number, k, of ancestors.
        // Include k=1, because this is a finite segment.
        for(int k=1; k <= n; ++k) {
            assert(0 == PtrLst_length(sd.a));
            sd.a = a[k-1];
            sd.nparts = k;
            sd.elen = elen[k-1];
            sd.lnconst = lnCoalConst(n, k);
            status = traverseSetPartitions(n, k, visitSetPart, &sd);
            if(status)
                return status;
        }
    }

    // max is maximum number of ancestors.
    int max;
    for(max = self->max; 0 == PtrLst_length(a[max-1]); --max) {
        fprintf(stderr,"%s:%d: deleting a[%d]\n",
                __FILE__,__LINE__, max-1);
        PtrLst_free(a[max-1]);
    }

    if(self->nparents == 1) {
        // Move all IdSet objects to parental waiting room.
        nw = self->parent[0]->nw;
        self->parent[0]->w[nw] = malloc(max * sizeof(PtrVec *));
        CHECKMEM(self->parent[0]->w[nw]);
        for(i=0; i<max; ++i) {
            self->parent[0]->w[nw] = PtrLst_to_PtrVec(a[i], NULL);
            PtrLst_free(a[i]);
        }
        self->parent[0]->nw += 1;
    }else {
        assert(self->nparents == 2);
        assert(*self->mig > 0.0);
        /*
          P[x] = (k choose x) * m^x * (1-m)^(k-x)
         */
        MigDat msd = {
                      .migrants = PtrLst_new(),
                      .natives = PtrLst_new(),
                      .a = NULL,
                      .migrationEvent = nextMigrationEvent();
        };
        for(k=1; k <= max; ++k) {
            msd.a = PtrLst_to_PtrVec(a[k-1], msd.a);
            PtrLst_free(a[k-1]);
            for(long x=0; x <= k; ++x) {
                // prob that x of k lineages are migrants
                long double lnpr = lbinomial(k, x);
                lnpr += x*logl(self->mix);
                lnpr += (k-x)*logl(1.0-self->mix);
                msd.pr = expl(lnpr);
                msd.nMigrants = x;
                msd.nNatives = k - x;
                status = traverseSetPartitions(n, k, visitMig, &msd);
                if(status)
                    return status;
            }
        }
    }

    // Free IdSet objects of descendants
    ids = PtrVec_pop(self->d[n-1]);
    while( ids ) {
        IdSet_free(ids);
        ids = PtrVec_pop(self->d[n-1]);
    }
    return status;
}

static int Segment_coalesceInfinite(Segment *self, double v, int dosing,
                                    BranchTab *branchtab) {
    double pr[self->max], elen[self->max];
    int n, i, status=0;
    IdSet *ids;

    CombDat cd = {.branchtab = branchtab,
                  .dosing = dosing
    };

    // Outer loop over numbers of descendants.
    // Calculate expected branch lengths, elen.
    for(n=2; n <= self->max; ++n) {

        // eigenvalues of transient states
        double eig[n-1];
        memset(eig, 0, (n-1)*sizeof(eig[0]));
        
        // Calculate expected length, within the segment, of
        // coalescent intervals with 2,3,...,n lineages.  elen[i] is
        // expected length of interval with i-1 lineages.
        // elen[0] not set.
        MatCoal_ciLen(n-1, elen+1, eig);

        // Calculate expected length, elen[0], of interval with one
        // lineage.  elen[i] refers to interval with i+1 lineages.
        double sum = 0.0;
        for(i=n; i >= 2; --i)
            sum += elen[i-1];
        elen[0] = v - sum; // elen[0] is infinite if v is.

        // Multiply elen by probability that segment had n lineages on
        // recent end of segment.
        for(i=0; i < n; ++i)
            elen[i] *= self->p[0][n-1];

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

        // Free IdSet objects of descendants
        ids = PtrVec_pop(self->d[n-1]);
        while( ids ) {
            IdSet_free(ids);
            ids = PtrVec_pop(self->d[n-1]);
        }
    }
    return status;
}

int Segment_coalesce(Segment *self, int dosing, BranchTab *branchtab) {
    double v;
    int n, i, status=0;
    IdSet *ids;

    if(self->nchildren > 0) {
        status = Segment_coalesce(self->child[0], dosing, branchtab);
        if(status)
            return status;
    }
    if(self->nchildren ==2 ) {
        status = Segment_coalesce(self->child[1], dosing, branchtab);
        if(status)
            return status;
    }

    if(self->end == NULL)
        v = INFINITY;
    else
        v = *self->end = *self->start;

    const int finite = isfinite(v); // is segment finite?

    if(finite)
        status = Segment_coalesceFinite(self, v, dosing, branchTab);
    else
        status = Segment_coalesceInfinite(self, v, dosing, branchTab);

    // Transfer ancestors to parent node or nodes
    if(self->nparents == 1) {
        for(i=0; i < self->max; ++i) {
            ids = PtrVec_pop(self->d[i]);
            while( ids ) {
                PtrVec_push(self->parent[0]->d[i], ids);
                ids = PtrVec_pop(self->d[i]);
            }
        }
    }else if(self->nparents == 2) {
    }
    
    return status;
}

/// Visit a combination
int visitComb(int d, int ndx[d], void *data) {
    assert(d>0);
    CombDat *dat = (CombDat *) data;

    unsigned nIdSets = PtrVec_length(self->d);
    for(unsigned i=0; i < nIdSets; ++i) {
        IdSet *ids = PtrVec_get(self->d, i);
        
        // sitepat is the union of the current set of descendants, as
        // described in ndx.
        tipId_t sitepat = 0;
        for(int i=0; i < d; ++i) {
            assert(ndx[i] < ids->nIds);
            sitepat |= ids->tid[ndx[i]];
        }
        
        // Skip singletons unless data->dosing is nonzero
        if(!dat->dosing && isPow2(sitepat))
            continue;
      
        // Increment BranchTab entry for current sitepat value.
        BranchTab_add(dat->branchtab, sitepat,
                      ids->p * dat->contrib);
    }
    return 0;
}

/// Visit a set partition.
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

    unsigned nIds = PtrVec_length(vdat->d);
    for(unsigned i=0; i < nIds; ++i) {
        IdSet *descendants = PtrVec_get(vdat->d, i);
        tipId_t sitepat[k];
        memset(sitepat, 0, k*sitepat[0]);

        assert(n == descendants->nIds);

        // Loop over descendants, creating a sitepat for each
        // ancestor. a[i] is the index of the ancestor of the i'th
        // descendant. sitepat[j] is the site pattern of the j'th
        // ancestor.
        for(int i=0; i<n; ++i)
            sitepat[a[i]] |= descendants->tid[i];

        // Loop over ancestors, i.e. over site patterns, adding
        // to the corresponding entry in BranchTab.
        for(int i=0; i<k; ++i)
            BranchTab_add(vdat->branchtab, sitepat[i],
                          p * vdat->elen * descendants->p);

        // Add the current set partition to the list of ancestral
        // states.
        IdSet *ancestors = IdSet_new(k, sitepat, p * descendants->p);
        IdSet_copyMigoutcome(ancestors, descendants);
        status = PtrLst_push(vdat->a, ancestors);
        if(status) {
            fprintf(stderr,"%s:%d can't push ancestors; status=%d\n",
                    __FILE__,__LINE__, status);
            exit(EXIT_FAILURE);
        }
    }

    return status;
}

/// Visit a set partition defining migrants
int visitMig(unsigned n, unsigned a[n], void *data) {
    MigDat *vdat = (MigDat *) data;
    int status=0;

    tipId_t migrants[vdat->nMigrants], natives[vdat->nNatives];

    // number of sets of ancestors
    int nSets = PtrVec_length(vdat->a);
    
    for(int i_set=0; i_set < nSets; ++i) {
        IdSet *set = PtrVec_get(vdat->a, i_set);

        int i_mig=0, i_nat=0;

        // tipId_t values of migrants and natives
        for(int i=0; i<n; ++i) {
            if(a[i]) // migrant
                migrants[i_mig++] = vdat->tid[i];
            else     // native
                natives[i_nat++] = vdat->tid[i];
        }

        // Perhaps this test should always be done
        assert(i_mig == vdat->nMigrants);
        assert(i_nat == vdat->nNatives);

        // Create IdSet objects for migrants and natives
        IdSet *mig = IdSet_new(i_mig, migrants, vdat->pr * set->p);
        IdSet *nat = IdSet_new(i_nat, natives, vdat->pr * set->p);
        PtrLst_push(vdat->migrants, mig)
        PtrLst_push(vdat->natives, nat)
    }

    return status;
}

/// If to==NULL, then return a new PtrVec, containing all the
/// pointers in "from". On return, "from" is empty, but the
/// PtrLst object itself is not freed.
///
/// If to!=NULL, then resize "to" so that it is large enough to hold
/// all the pointers in "from", and move all pointers from "from" into
/// "to". Return "to".
PtrVec *PtrLst_to_PtrVec(PtrLst *from, PtrVec *to) {
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
    
    IdSet *idset = PtrLst_pop(from);
    while( idset ) {
        PtrVec_push(to, idset);
        idset = PtrLst_pop(from);
    }

    return to;
}
