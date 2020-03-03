#include "segment.h"
#include "partprob.h"
#include "intpart.h"
#include "comb.h"
#include "binary.h"
#include "matcoal.h"
#include "misc.h"
#include "branchtab.h"
#include <stdlib.h>

typedef struct IdSet IdSet;
typedef struct CombDat CombDat;

struct CombDat {
    double contribution; // Pr[site pattern]*E[len of interval]
    IdSet *ids;
    BranchTab *branchtab;
    int dosing;    // do singleton site patterns if nonzero
};

// A set of tipId_t values.
struct IdSet {
    int nIds; // number of Ids in this set
    double p; // probability of this set

    // Array of length nIds. tid[i] is the id of the i'th haploid
    // individual in this set.
    tipId_t *tid;

    // A tipIt_t value representing the union of all individuals in
    // this set. It is illegal to join two IdSet objects, a and b, if
    // they overlap. In other words, if "a.allbits & b.allbits" is not
    // zero.  The i'th bit is set in allbits if it is set in any of
    // the constituent tipId_t values.
    tipId_t allbits;

    IdSet *next;
};

// Probabilistic variables for one segment of a population tree, i.e.,
// for ont PopNode.
struct Segment {
    int allocated; // allocated dimension of arrays
    int max;     // max number of lineages in segment

    // p[0][i] is prob there are i+1 lineages at recent end of segment
    // p[1][i] is analogous prob for ancient end of interval.
    double *p[2];

    // Arrays of pointers to linked lists of IdSet objects. Dimension
    // is max X 2.  ids[0] refers to the recent end of the segment and
    // ids[1] to the ancient end. ids[0][i] is the list for the case
    // in which there are i+1 lineages at the recent end of the segment.
    IdSet **ids[2];
};

int tipId_vec_cmp(int nx, const tipId_t x[nx], int ny,
                  const tipId_t y[ny]);
IdSet *IdSet_new(IdSet *next, int nIds, tipId_t tid[nIds], double prob);
IdSet *IdSet_add(IdSet *head, int nIds, tipId_t tid[nIds], double prob);
void IdSet_free(IdSet *self);
IdSet *IdSet_join(IdSet *a, IdSet *b);
int IdSet_length(IdSet *self);
void IdSet_mulBy(IdSet *self, double factor);
void IdSet_divBy(IdSet *self, double divisor);
double IdSet_sumProb(IdSet *self);
int visitComb(int d, int ndx[d], void *data);

/**
 * Compare two vectors, x and y, of tipId_t values. Return -1 if x<y,
 * 1 if x>y, and 0 otherwise. 
 *
 * @param[in] nx length of x
 * @param[in] x  vector of tipId_t values
 * @param[in] ny length of y
 * @param[in] y vector of tipId_t values
 */
int tipId_vec_cmp(int nx, const tipId_t x[nx], int ny,
                   const tipId_t y[ny]) {
    for(int i=0; i < nx && i < ny; ++i) {
        if(x[i] < y[i])
            return -1;
        if(x[i] > y[i])
            return 1;
    }
    if(nx < ny)
        return -1;
    if(nx > ny)
        return 1;
    return 0;
}

/**
 * Allocate a new IdSet object with given values of nIds, tid, and
 * prob, and with the next field pointing to the argument "next".
 */
IdSet *IdSet_new(IdSet *next, int nIds, tipId_t tid[nIds], double prob) {
    IdSet *self = malloc(sizeof(IdSet));
    CHECKMEM(self);

    self->nIds = nIds;
    self->p = prob;
    self->tid = malloc(nIds * sizeof(self->tid[0]));
    CHECKMEM(self->tid);

    self->allbits = 0;
    for(int i=0; i < nIds; ++i) {
        self->allbits |= tid[i];
        self->tid[i] = tid[i];
    }

    self->next = next;
    return self;
}

/**
 * Add a new IdSet item to a sorted list. If a corresponding IdSet
 * object already exists, then add prob to that object. The value
 * of nIds must match the value in the existing list.
 */
IdSet *IdSet_add(IdSet *head, int nIds, tipId_t tid[nIds], double prob) {
    if(head==NULL)
        return IdSet_new(NULL, nIds, tid, prob);
    if(head->nIds != nIds) {
        fprint(stderr,"%s:%d: can't add a set with %d ids to a list"
               " of sets with %d ids\n",
               __FILE__,__LINE__, nIds, head->nIds);
        exit(EXIT_FAILURE);
    }
    int cmp = tipId_vec_cmp(nIds, tid, head->nIds, head->tid);
    if(cmp < 0)
        return IdSet_new(head, nIds, tid, prob);
    if(cmp > 0) {
        head->next = IdSet_add(head->next, nIds, tid, prob);
        return head;
    }
    head->p += prob;
    return head;
}

void IdSet_free(IdSet *self) {
    if(self==NULL)
        return;
    IdSet_free(self->next);
    free(self->tid);
    free(self);
}

/**
 * Join two sets. The probability of the joined set is the product
 * of the two input probabilities. The joined set is the union of
 * sets a and b. If sets a and b overlap, the join is not possible,
 * and the function returns NULL. This happens when we are trying to
 * join sets that represent mutually exclusive outcomes. For example,
 * as we trace a nucleotide's history backwards through the network of
 * populations, we come to branch points. One branch represents
 * migration from another population, the other represents descent
 * from the same population. These are mutually exclusive outcomes,
 * which both occur with nonzero probability. Eventually, as we move
 * farther back in time, the algorithm will try to join the ancestors
 * of these two outcomes. At that point, IdSet_join will return NULL,
 * and we'll know that this pair of sets cannot be joined.
 */
IdSet *IdSet_join(IdSet *a, IdSet *b) {
    if( (a->allbits & b->allbits) != 0) {
        // Sets overlap. Don't join them.
        return NULL;
    }

    IdSet *self = malloc(sizeof(IdSet));
    CHECKMEM(self);

    // sum
    self->nIds = a->nIds + b->nIds;

    // product
    self->p = a->p * b->p;

    // The i'th bit is set in allbits if it is set in any of
    // the constituent tipId_t values.
    self->allbits = a->allbits | b->allbits;

    // New array is large enough to hold both sets of tipId_t values.
    self->tid = malloc(self->nIds * sizeof(self->tid[0]));
    CHECKMEM(self->tid);

    // Copy arrays into self, maintaining sort.
    int ia=0, ib=0, j=0;
    while(ia < a->nIds && ib < b->nIds) {
        assert(j < self->nIds);
        if(a->tid[ia] < b->tid[ib])
            self->tid[j++] = a->tid[ia++];
        else {
            assert(a->tid[ia] > b->tid[ib]);
            // a->tid can't equal b->tid, because the first step in
            // this function ensures that no bit is set in both
            // values. The assertion just above checks this.
            self->tid[j++] = b->tid[ib++];
        }
    }
    while(ia < a->nIds)
        self->tid[j++] = a->tid[ia++];
    while(ib < b->nIds)
        self->tid[j++] = a->tid[ib++];
    assert(ia == a->nIds);
    assert(ib == b->nIds);
    assert(j == self->nIds);

    return self;
}

int IdSet_length(IdSet *self) {
    int len=0;
    while(self != NULL) {
        len += 1;
        self = self->next;
    }
    return len;
}

void IdSet_mulBy(IdSet *self, double factor) {
    while(self != NULL) {
        self->p *= factor;
        self = self->next;
    }
}

void IdSet_divBy(IdSet *self, double divisor) {
    while(self != NULL) {
        self->p /= divisor;
        self = self->next;
    }
}

double IdSet_sumProb(IdSet *self) {
    double sum=0.0;
    while(self != NULL) {
        sum += self->p;
        self = self->next;
    }
    return sum;
}

Segment *Segment_new(void) {
    Segment *self = malloc(sizeof(Segment));
    CHECKMEM(self);

    self->allocated = 8;
    self->max = 0;

    self->p[0] = malloc(self->allocated * sizeof(self->p[0][0]));
    CHECKMEM(self->p[0]);

    self->p[1] = malloc(self->allocated * sizeof(self->p[1][0]));
    CHECKMEM(self->p[1]);

    self->ids[0] = self->ids[1] = NULL;

    return self;
}

// Add IdSet to Segment
void Segment_add(Segment *self, IdSet *idset) {
    int nids = idset->nids;

    // reallocate arrays if necessary
    if(nids > self->allocated) {
        if(2*self->allocated >= nids)
            self->allocated *= 2;
        else
            self->allocated = nids;

        self->p[0] = realloc(self->p[0], self->allocated * sizeof(self->p[0][0]));
        CHECKMEM(self->p[0]);
        self->p[1] = realloc(self->p[1], self->allocated * sizeof(self->p[1][0]));
        CHECKMEM(self->p[1]);

        self->ids[0] = realloc(self->ids[0], self->allocated * sizeof(self->ids[0][0]));
        CHECKMEM(self->ids[0]);
        self->ids[1] = realloc(self->ids[1], self->allocated * sizeof(self->ids[1][0]));
        CHECKMEM(self->ids[1]);
    }

    // add IdSet to segment
    self->ids[0][nids-1] = IdSet_add(self->ids[0][nids-1], nids, idset->tid, 0.0);
}

int Segment_coalesce(Segment *self, int maxsamp, int dosing,
                     BranchTab *branchtab, double v) {
    assert(self->max <= maxsamp);

    double pr[maxsamp], elen[maxsamp], sum;
    int n, i, status=0;
    const int finite = isfinite(v); // is segment finite?
    int kstart = (finite ? 1 : 2);  // skip k=1 if infinite

    // If there is only one line of descent, no coalescent events are
    // possible, so p[1][0] is at least as large as p[0][0].
    self->p[1][0] = self->p[0][0];

    // Calculate probabilities and expected values, p[1] and elen.
    for(n=2; n <= self->max; ++n) {

        // Skip improbable states.
        if(self->p[0][n-1] == 0.0)
            continue;
        
        // Calculate prob of 2,3,...,n lines of descent.
        // On return, pr[1] = prob[2], pr[n-1] = prob[n],
        // pr[0] not set.
        MatCoal_project(n-1, pr+1, v);

        // Calculate probability of 1 line of descent. I'm doing the
        // sum in reverse order on the assumption that higher indices
        // will often have smaller probabilities.
        // pr[i] is prob of i+1 lineages; pr[0] not set
        sum=0.0;
        for(i=n; i > 1; --i)
            sum += pr[i-1];
        self->p[1][0] += self->p[0][n-1] * (1.0-sum);

        // Add probs of 2..n lines of descent
        // pr[i] is prob of i+2 lineages
        // self->p[1][i] is prob if i+1 lineages
        for(i=2; i <= n; ++i)
            self->p[1][i-1] += self->p[0][n-1] * pr[i-1];

        // Calculate expected length, within the segment, of
        // coalescent intervals with 2,3,...,n lineages.  elen[i] is
        // expected length of interval with i-1 lineages.
        // elen[0] not set.
        MatCoal_ciLen(n-1, elen+1, v);

        // Calculate expected length, elen[0], of interval with one
        // lineage.  elen[i] refers to interval with i+1 lineages.
        sum = 0.0;
        for(i=n; i >= 2; --i)
            sum += elen[i-1];
        elen[0] = v - sum; // elen[0] is infinite if v is.

        // Multiply elen by probability that segment had n lineages on
        // recent end of segment.
        for(i=1; i <= n; ++i)
            elen[i-1] *= self->p[0][n-1];
    }

    CombDat cd = {.contribution = 0.0,
                  .ids = NULL,
                  .branchtab = branchtab,
                  .dosing = dosing
    };

    // loop over number of descendants in this segment
    for(n=1; n <= self->max; ++n) {
        // Skip improbable states.
        if(self->p[0][n-1] == 0.0)
            continue;

        cd.ids = ids[0][n-1];

        // loop over intervals: k is the number of ancestors
        // w/i interval.
        for(int k=kstart; k <= n; ++k) {
            
            // portion of log Qdk that doesn't involve d
            long double lnconst = logl(k) - lbinom(n-1, k-1);

            // Within each interval, there can be ancestors
            // with 1 descendant, 2, 3, ..., n-k+2.
            for(int d=1; d <= n-k+1; ++d) {
                long double lnprob = lnconst
                    + lbinom(n-d-1, k-1) - lbinom(n,d);

                // probability of site pattern
                cd.contribution = (double) expl(lnprob);

                // times expected length of interval
                cd.contribution *= elen[k-1];
                
                status = traverseComb(n, d, visitComb, &cd);
                if(status)
                    return status;
            }
        }
    }
    return status;
}

/// Visit a combination
int visitComb(int d, int ndx[d], void *data) {
    assert(d>0);
    CombDat *dat = (CombDat *) data;

    for(IdSet ids = dat->ids; ids != NULL; ids = ids->next) {
        
        // tid is the union of descendants in current set.
        tipId_t tid = 0;
        for(int i=0; i < d; ++i)
            tid |= ids->tid[ndx[i]];

        // Skip singletons unless data->dosing is nonzero
        if(!dat->dosing && isPow2(tid[i]))
            continue;
      
        // Increment BranchTab entry for current tid value.
        BranchTab_add(dat->branchtab, tid,
                      ids->p * dat->contribution);
    }
    return 0;
}
