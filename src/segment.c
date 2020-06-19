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
#include "nodestore.h"
#include "segment.h"
#include <stdlib.h>
#include <string.h>

typedef struct MigDat MigDat;
typedef struct CombDat CombDat;
typedef struct SetPartDat SetPartDat;

// For assigning initial tipId_t values within segments.
static tipIt_t currTipId = 0;

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
    PtrVec *a;

};

// Data manipulated by visitSetPart function.
struct SetPartDat {
    unsigned nparts; // Number of parts in partition
    long double lnconst;  // log of constant in Durrett's theorem 1.5
    double elen;     // E[len of interval]

    // PtrVec_get(a, i) is the IdSet for the i't set of ancestors
    PtrVec *a;

    // PtrVec_get(d, i) is the IdSet for the i'th set of descendants
    PtrVec *d;

    BranchTab *branchtab;
    int dosing;      // do singleton site patterns if nonzero
};

void   Segment_addIdSet(Segment *self, IdSet *idset);
int    visitComb(int d, int ndx[d], void *data);
int    visitSetPart(unsigned n, unsigned a[n], void *data);

// The total number of samples.
static int total_samples = 0;

void *Segment_new(double *twoN, double *start, int nsamples,
                  NodeStore * ns) {
    total_samples += nsamples;
    
    Segment *self = NodeStore_alloc(ns);
    CHECKMEM(self);

    memset(self, 0, sizeof(*self));

    self->twoN = twoN;
    self->start = start;

    self->nsamples = nsamples;
    for(int i=0; i < nsamples; ++i)
        self->sample[i] = currTipId++;

    // to be set later: max, a, d, nw, wdim, w.

    return self;
}

int      Segment_addChild(void * vparent, void * vchild) {
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
    if(*child->start > *parent->start) {
        fprintf(stderr,
                "%s:%s:%d: Child start (%lf) must be <= parent start (%lf)\n",
                __FILE__, __func__, __LINE__, *child->start, *parent->start);
        return DATE_MISMATCH;
    }
    if(child->end == NULL) {
        child->end = parent->start;
    } else {
        if(child->end != parent->start) {
            fprintf(stderr, "%s:%s:%d: Date mismatch."
                    " child->end=%p != %p = parent->start\n",
                    __FILE__, __func__, __LINE__, child->end, parent->start);
            return DATE_MISMATCH;
        }
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
    REQUIRE(self->twoN != NULL, file, lineno);
    REQUIRE(*self->twoN > 0.0, file, lineno);
    REQUIRE(self->start != NULL, file, lineno);
    REQUIRE(*self->start >= 0.0, file, lineno);
    if(self->end) {
        REQUIRE(self->nparents > 0, file, lineno);
        REQUIRE(*self->start <= *self->end, file, lineno);
    }
    switch(self->nparents) {
    case 0:
        REQUIRE(self->end == NULL, file, lineno);
        REQUIRE(self->mix == NULL, file, lineno);
        break;
    case 1:
        REQUIRE(self->end != NULL, file, lineno);
        REQUIRE(self->mix == NULL, file, lineno);
        break;
    case 2:
        REQUIRE(self->end != NULL, file, lineno);
        REQUIRE(self->mix != NULL, file, lineno);
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
    if(self->nchildren == 2)
        Segment_sanityCheck(self->child[1], file, lineno);
#endif
}

int      Segment_mix(void * vchild, double *mPtr, void * vintrogressor, 
                     void * vnative) {
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
    if(child->end != NULL) {
        if(child->end != introgressor->start) {
            fprintf(stderr,"%s:%s:%d: Date mismatch."
                    " child->end=%p != %p=introgressor->start\n",
                    __FILE__, __func__, __LINE__,
                    child->end, introgressor->start);
            return DATE_MISMATCH;
        }
        if(child->end != native->start) {
            fprintf(stderr, "%s:%s:%d: Date mismatch."
                    " child->end=%p != %p=native->start\n",
                    __FILE__, __func__, __LINE__, child->end, native->start);
            return DATE_MISMATCH;
        }
    } else if(native->start != introgressor->start) {
        fprintf(stderr, "%s:%s:%d: Date mismatch."
                "native->start=%p != %p=introgressor->start\n",
                __FILE__, __func__, __LINE__,
                native->start, introgressor->start);
        return DATE_MISMATCH;
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
    Segment_sanityCheck(child, __FILE__, __LINE__);
    Segment_sanityCheck(introgressor, __FILE__, __LINE__);
    Segment_sanityCheck(native, __FILE__, __LINE__);
    return 0;
}

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
    fprintf(fp, "%p twoN=%lf ntrval=(%lf,", self, *self->twoN, *self->start);
    if(self->end != NULL)
        fprintf(fp, "%lf)\n", *self->end);
    else
        fprintf(fp, "Inf)\n");

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

int Segment_coalesce(Segment *self, int dosing, BranchTab *branchtab) {
    double v, pr[self->max], elen[self->max];
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

    CombDat cd = {.branchtab = branchtab,
                  .dosing = dosing
    };

    SetPartDat sd = {.branchtab = branchtab,
                     .dosing = dosing
    };

#ifndef NDEBUG
    for(i=0; i < self->max; ++i)
        assert(self->ids[1][i] == NULL);
#endif    

    // We only need to calculate the probabilities, self->p[1][i],
    // of ancestors if the segment is finite.
    if(finite) {
        memset(self->p[1], 0, MAXSAMP*sizeof(double));

        // If there is only one line of descent, no coalescent events
        // are possible, so p[1][0] is at least as large as p[0][0].
        self->p[1][0] = self->p[0][0];
    }

    // Outer loop over numbers of descendants.
    // Calculate probabilities and expected values, p[1] and elen.
    for(n=2; n <= self->max; ++n) {

        // eigenvalues of transient states
        double eig[n-1];
        if(finite)
            MatCoal_eigenvals(n-1, eig, v);
        else
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

        if(finite) {
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

            sd.a = self->a[n-1];

            // Loop over number, k, of ancestors.
            // Include k=1, because this is a finite segment.
            for(int k=1; k <= n; ++k) {
                sd.d = self->d[k-1];
                sd.nparts = k;
                sd.elen = elen[k-1];
                sd.lnconst = lnCoalConst(n, k);
                status = traverseSetPartitions(n, k, visitSetPart, &sd);
                if(status)
                    return status;
            }
            sd.a = sd.d = NULL;
        }else{  // Infinite segment: the root

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

        // Free IdSet objects of descendants
        ids = PtrVec_pop(self->d[n-1]);
        while( ids ) {
            IdSet_free(ids);
            ids = PtrVec_pop(self->d[n-1]);
        }
    }

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
        status = PtrVec_push(vdat->a, ancestors);
        if(status) {
            fprintf(stderr,"%s:%d can't push ancestors; status=%d\n",
                    __FILE__,__LINE__, status);
            exit(EXIT_FAILURE);
        }
    }

    return status;
}
