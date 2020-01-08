#include "segment.h"
#include "intpart.h"
#include "comb.h"
#include "binary.h"
#include "matcoal.h"
#include "misc.h"
#include "branchtab.h"
#include <stdlib.h>

typedef struct IdSet IdSet;
typedef struct PartDat PartDat;

struct PartDat {
    int nsub; // number of subsets in partition
    double lnconst; // log of constant in Durrett's theorem 1.5
    double intLen; // expected length of this coalescent interval
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
    int max;     // max number of lineages in segment

    // p[0][i] is prob there are i+1 lineages at recent end of segment
    // p[1][i] is analogous prob for ancient end of interval.
    double *p[2];

    // Arrays of pointers to linked lists of IdSet objects. Dimension
    // is max X 2.  ids[0] refers to the recent end of the segment and
    // ids[1] to the ancient end. ids[0][i] is the list for the case
    // in which there are i+1 lineages.
    IdSet **ids[2];
};

long double lnCoalConst(unsigned n, unsigned k);
double probPartition(unsigned k, unsigned y[k], long double lnconst);
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
int visitPart(int k, int y[k], void *data);

/** 
 * Log of constant in coalescent probability from theorem 1.5, p. 11, Durrett,
 * Richard. 2008. Probability Models for DNA Sequence Evolution.
 */
long double lnCoalConst(unsigned n, unsigned k) {
    assert(n >= k);
    assert(k > 0);
    return lgammal(k+1)
        - lgammal(n+1)
        + lgammal(n-k+1)
        + lgammal(k)
        - lgammal(n);
}

/**
 * Partition probability under the coalescent process. There are n
 * descendants in some recent epoch and k < n in some earlier
 * epoch. In that earlier epoch, the i'th ancestor had y[i]
 * descendants. The function returns the probability of a partition
 * that satisfies this condition, given n and k. There may be several
 * such partitions. lnconst should be calculated using function
 * lnCoalConst. See theorem 1.5, p. 11, Durrett,
 * Richard. 2008. Probability Models for DNA Sequence Evolution.
 */
double probPartition(unsigned k, unsigned y[k], double lnconst) {
    assert(k > 0);
    long double x = 0.0L;
    for(unsigned i=0; i<k; ++i) {
        x += lgammal(y[i] + 1);
    }
    return (double) expl(lnconst + x);
}

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
 * object already exists, then add prob to that object.
 */
IdSet *IdSet_add(IdSet *head, int nIds, tipId_t tid[nIds], double prob) {
    if(head==NULL)
        return IdSet_new(NULL, nIds, tid, prob);
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

int Segment_coalesce(Segment *self, int maxsamp, int dosing,
                     BranchTab *branchtab, double v) {
    assert(self->max <= maxsamp);

    double x[maxsamp], sum;
    int n, i, status=0;
    
    // If there is only one line of descent, no coalescent events are
    // possible, so p[1][0] is at least as large as p[0][0].
    self->p[1][0] = self->p[0][0];

    // Calculate probabilities and expected values, p[1] and cilen.
    for(n=2; n <= self->max; ++n) {

        // Skip improbable states.
        if(self->p[0][n-1] == 0.0)
            continue;
        
        // Calculate prob of 2,3,...,n lines of descent.
        // On return, x[1] = prob[2], x[n-1] = prob[n],
        // x[0] not set.
        MatCoal_project(n-1, x+1, v);

        // Calculate probability of 1 line of descent. I'm doing the
        // sum in reverse order on the assumption that higher indices
        // will often have smaller probabilities.
        // x[i] is prob of i+1 lineages; x[0] not set
        sum=0.0;
        for(i=n; i >= 2; --i)
            sum += x[i-1];
        self->p[1][0] += self->p[0][n-1] * (1.0-sum);

        // Add probs of 2..n lines of descent
        // x[i] is prob of i+2 lineages
        // self->p[1][i] is prob if i+1 lineages
        for(i=2; i <= n; ++i)
            self->p[1][i-1] += self->p[0][n-1] * x[i-1];

        // Calculate expected length, within the segment, of
        // coalescent intervals with 2,3,...,n lineages.  x[i] is
        // expected length of interval with i-1 lineages.
        // x[0] not set.
        MatCoal_ciLen(n-1, x+1, v);

        // Calculate expected length, x[0], of interval with one
        // lineage.  x[i] refers to interval with i+1 lineages.
        sum = 0.0;
        for(i=n; i >= 2; --i)
            sum += x[i-1];
        x[0] = v - sum;

        // Multiply x by probability that segment had n lineages on
        // recent end of segment.
        for(i=1; i <= n; ++i)
            x[i-1] *= self->p[0][n-1];

        for(int k=1; k <= n; ++k) {
            PartDat pd = {.lnconst = lnCoalConst(n, k),
                          .nsub = k,
                          .intLen=x[k-1],
                          .ids = self->ids[0][k-1],
                          .branchtab = branchtab,
                          .dosing = dosing};
            status = traverseSetPartitions(n, k, visitPart, &pd);
            if(status)
                return status;
        }
    }
    return status;
}

/// Visit a set partition.
int visitPart(int n, int y[n], void *data) {
    PartDat *dat = (PartDat *) data;

    // c[i] is number of descendants in subset i
    int c[dat->nsub];
    memset(c, 0, dat->nsub * sizeof(c[0]));
    for(int i=0; i < n; ++i) {
        assert(y[i] < dat->nsub);
        ++c[y[i]];
    }

    // probability inherited from upstream
    double x = dat->prob;

    // times probability of this partition
    x *= probPartition(dat->nsub, c, dat->lnconst);

    // times expected length of interval
    x *= dat->intLen;

    // tid[i] is the union of the descendants of ancestor i.
    tipId_t tid[dat->nsub];
    memset(tid, 0, dat->nsub * sizeof(tid[0]));
    for(int i=0; i < n; ++i)
        tid[y[i]] |= dat->tid[i];

    // For each ancestor, increment the corresponding
    // branchtab entry.
    for(int i=0; i < dat->nsub; ++i) {
        // Skip singletons unless data->dosing is nonzero
        if(!dat->dosing && isPow2(tid[i]))
            continue;
      
        // Increment BranchTab entry for current tid value.
        BranchTab_add(dat->branchtab, tid, x);
    }
    return 0;
}
