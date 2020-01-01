#include "segment.h"
#include "intpart.h"
#include "comb.h"
#include "binary.h"

typedef struct IdSet IdSet;
typedef struct IdSetList IdSetList;
typedef struct PartDat PartDat;

struct PartDat {
    double prob;   // reciprocal of (n-1) choose (k-1)
    double prcomb; // prob of each allocation w/i current partition
    double intLen; // expected length of this coalescent interval
    IdSet *head;
    BranchTab *bt;
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
    tipIt_t allbits;

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
int visitIntPart(int k, int y[k], void *data);
int visitComb(int k, int y[k], int *ndx[k], void *data);

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
    head->prob += prob;
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

int Segment_coalesce(Segment *self, int maxsamp,
                     MatCoal *mc[maxsamp-1],
                     Stirling2 *stirling2,
                     double v) {
    assert(self->max <= maxsamp);

    double x[maxsamp-1], sum;
    int n, i;
    
    // If there is only one line of descent, no coalescent events are
    // possible, so p[1][0] is at least as large as p[0][0].
    self->p[1][0] = self->p[0][0];

    // Calculate probabilities and expected values, p[1] and cilen.
    for(n=2; n <= self->max; ++n) {

        // Skip improbable states.
        if(self->p[0][n-1] == 0.0)
            continue;
        
        // Calculate prob of 2,3,...,n lines of descent.
        // On return, x[0] = prob[2], x[n-2] = prob[n]
        MatCoal_project(n-1, x, v);

        // Calculate probability of 1 line of descent. I'm doing the
        // sum in reverse order on the assumption that higher indices
        // will often have smaller probabilities.
        // x[i] is prob of i+2 lineages
        sum=0.0;
        for(i=n; i >= 2; --i)
            sum += x[i-2];
        self->p[1][0] += self->p[0][n-1] * (1.0-sum);

        // Add probs of 2..n lines of descent
        // x[i] is prob of i+2 lineages
        // self->p[1][i] is prob if i+1 lineages
        for(i=2; i <= n; ++i)
            self->p[1][i-1] += self->p[0][n-1] * x[i-2];

        // Calculate expected length, within the segment, of
        // coalescent intervals with 2,3,...,n lineages.  x[i] is
        // expected length of interval with i-2 lineages.
        MatCoal_ciLen(n-1, x, v);

        // Calculate expected length of interval with one lineage.
        // x[i] refers to interval with i+2 lineages.
        sum = 0.0;
        for(i=n; i >= 2; --i)
            sum += x[i-2];
        double cilen = self->pr[0][n-1] * (v - sum);

        // Add lengths of intervals with 2,3,...,n lineages.
        // x[i] refers to interval with i+2 lineages.
        // self->ci[1][i] refers to intervals with i+1 lineages.
        for(i=2; i <= n; ++i)
            self->ci[i-1] += self->pr[0][n-1] * x[i-2];
    }
}

/// Visit an integer partition.
int visitIntPart(int k, int y[k], void *data) {
    PartDat *dat = (PartDat *) data;

    int i, j;

    // c[0] is number of y[i] with largest value, c[1] is number
    // with next largest value, etc. Algorithm relies of fact
    // that entries of y are sorted.
    int c[k];
    int m=0; // number of distinct values in partition y
    c[0] = 1;
    for(j=1; j<k; ++j) {
        if(y[j] == y[j-1])
            ++c[m];
        else
            c[++m] = 1;
    }

    // Number of x vectors contributing to this partition.
    long n = multinom(m, c);
    
    dat->prcomb = nx * dat->prob;  // prob of current partition

    n = multinom(k, y);
    dat->prcomb /= n;              // prob of each combination

    // Visit all ways of allocating descendants to ancestors.
    return traverseMultiComb(k, y[k], visitComb, data);
}

/**
 * Visit a particular allocation of descendants among ancestors, and
 * increment the BranchTab object.
 * @param[in] k The number of ancestors.
 * @param[in] y y[i] is the number of descendants of ancestor i. The
 * sum of y equals n, the total number of descendants.
 * @param[in] ndx ndx[i][j] is the index of the j'th descendant of
 * ancestor i. 
 * Each index is between 0 and n-1 and appears only once in the
 * two-dimensional array, ndx.
 * @param[inout] data Pointer to an object of type PartDat.
 */
int visitComb(int k, int y[k], int *ndx[k], void *data) {
    PartDat *dat = (PartDat *) data;

    // Loop across sets of ancestors. These sets may have different
    // probabilities, as indicated by s->p.
    for(IdSet *s = dat->ids; s!=NULL; s = s->next) {

        // Loop across ancestors
        for(int i=0; i<k; ++i) {

            // Loop across descendants of ancestor i to create a
            // tipId_t value representing the union of that ancestor's
            // descendants.
            tipId_t tid = 0;
            for(int j=0; j < y[i]; ++j) {
                assert(ndx[j] < s->nIds);
                tid |= s->tid[ndx[i][j]];
            }

            assert(tid);

            // Skip singletons unless data->dosing is nonzero
            if(!dat->dosing && isPow2(tid))
                continue;

            // probability times expected length of interval
            double x = s->p * dat->prcomb * dat->intLen;

            // Increment BranchTab entry for current tid value.
            BranchTab_add(dat->bt, tid, x);
        }
    }
    return 0;
}
