#ifndef ARR_IDSET_H
#  define ARR_IDSET_H

#  include "typedefs.h"
#  include <stdio.h>

// A set of tipId_t values.
struct IdSet {
    double p; // probability of this set
    int nIds; // number of Ids in this set

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

IdSet *IdSet_newTip(tipId_t tid);
IdSet *IdSet_new(IdSet *next, int nIds, tipId_t tid[nIds], double prob);
int    IdSet_cmp(const IdSet *x, const IdSet *y);
int    IdSet_deepCmp(const IdSet *x, IdSet *y);
IdSet *IdSet_dup(const IdSet *old);
IdSet *IdSet_add(IdSet *head, const IdSet *to_add, double prob);
void   IdSet_free(IdSet *self);
IdSet *IdSet_join(IdSet *a, IdSet *b);
int    IdSet_nIds(IdSet *self);
int    IdSet_nSets(IdSet *self);
int    IdSet_nDescendants(IdSet *self);
void   IdSet_mulBy(IdSet *self, double factor);
void   IdSet_divBy(IdSet *self, double divisor);
double IdSet_sumProb(IdSet *self);
void   IdSet_print(IdSet *self, FILE *fp);
void   IdSet_sanityCheck(IdSet *self, const char *file, int lineno);

#endif
