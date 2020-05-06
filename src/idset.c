/**
 * @file idset.c
 * @author Alan R. Rogers
 * @brief A set of tipID_t values.
 *
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "idset.h"

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

/**
 * Compare IdSet objects. Return -1 if x<y, 1 if x>y, and 0 otherwise.
 */
int IdSet_cmp(const IdSet *x, const IdSet *y) {
    for(int i=0; i < x->nIds && i < y->nIds; ++i) {
        if(x->tid[i] < y->tid[i])
            return -1;
        if(x->tid[i] > y->tid[i])
            return 1;
    }
    if(x->nIds < y->nIds)
        return -1;
    if(x->nIds > y->nIds)
        return 1;
    return 0;
}

/// Allocate a new IdSet with a single tipId_t value and probability 1.
IdSet *IdSet_newTip(tipId_t tid) {
    IdSet *self = malloc(sizeof(IdSet));
    CHECKMEM(self);

    self->nIds = 1;
    self->p = 1.0;
    self->tid = malloc(sizeof(self->tid[0]));
    CHECKMEM(self->tid);

    self->allbits = self->tid[0] = tid;

    self->next = NULL;
    return self;
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

IdSet *IdSet_dup(IdSet *old) {
    return IdSet_new(NULL, old->nIds, old->tid, old->prob);
}

/**
 * Add a new IdSet item to a sorted list. If a corresponding IdSet
 * object already exists, then add prob to that object. The value
 * of nIds must match the value in the existing list.
 */
IdSet *IdSet_add(IdSet *head, const IdSet *to_add, double prob) {
    IdSet *tmp;
    if(head==NULL) {
        tmp = IdSet_dup(to_add);
        tmp->prob *= prob;
        return tmp;
    }
    if(head->nIds != to_add->nIds) {
        fprint(stderr,"%s:%d: can't add a set with %d ids to a list"
               " of sets with %d ids\n",
               __FILE__,__LINE__, to_add->nIds, head->nIds);
        exit(EXIT_FAILURE);
    }
    int cmp = IdSet_cmp(to_add, head->nIds);
    if(cmp < 0) {
        tmp = IdSet_dup(to_add);
        tmp->prob *= prob;
        tmp->next = head;
        return tmp;
    }if(cmp > 0) {
        head->next = IdSet_add(head->next, to_add, prob);
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
