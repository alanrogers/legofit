#ifndef ARR_MIGOUTCOME
#define ARR_MIGOUTCOME

#include "typedefs.h"
#include <stdio.h>
#include <stdint.h>

/**
This structure handles bookkeeping associated with migration. A when
migration occurs, a number of mutually exclusive outcomes are
possible, and only ancestors of the same outcome can be allowed to
coalesce farther back in time.

Let k represent the number of ancestors at the ancient end of a
segment. Each of these is an immigrant with some positive
probability. There are 2^k mutually-exclusive ways to choose
immigrants and natives, and each of these outcomes divides the
ancestors into two sets, one of which may be empty.

In this structure, "event" is an index that indicates which episode of
migration is being described. It is 0 for the first migration event, 1
for the next, and so on.

"Outcome" is an integer whose i'th bit is on if this is the i'th
outcome.  It is forbidden to coalesce lineages that belong to
different outcomes of the same event.

When two populations coalesce, lineages in one may have a history of
migration events that didn't occur in the other. In that case, the
lineages of the lineages lacking this history are duplicated. One of
the two duplicate copies is associated with the lineages whose history
includes migration. The others acquire a new MigOutcome object whose
"event" and "noutcomes" field match those of the others. The "outcome"
field contains the complement of the "outcome" field of the other
lineages.

This complement can be calculated as

  y = ~0LU >> (8*sizeof(y) - noutcomes);
  y ^= x;

Here, the 1st line turns on all bits from 0 through noutcomes-1,
and the 2nd line turns off the bits that are on in x.

This results in MigOutcome objects in which multiple bits are on. In
effect, the corresponding lineages have been duplicated and are
associated with multiple outcomes of the migration event. They are not
really duplicated, however, until a merger of two populations (in
backwards time) brings them into contact with lineages with which they
share bits in the "outcome" field. At this point, the lineages are
duplicated. One copy acquires the bit pattern of the lineages with
which they are combined, and the other copy has these bits turned off.

A single set of lineages may have experienced multiple migration
events, so each lineages has a sorted linked list of MigOutcome
objects, which represent all the migration events it has experienced.
 **/
struct MigOutcome {
    unsigned event;        // which episode of migration?
    unsigned noutcomes;    // number of outcomes of this migration
                           // event

    // Bit i is on if this is the i'th outcome.
    uint64_t bits;

    double pr;             // probability of this outcome
    struct MigOutcome *next;
};

unsigned    nextMigrationEvent(void);
MigOutcome *MigOutcome_insert(MigOutcome *head,
                              unsigned event,
                              unsigned noutcomes,  
                              unsigned outcome,
                              double pr);
MigOutcome *MigOutcome_dup(MigOutcome *old);
void        MigOutcome_free(MigOutcome *self);
void        MigOutcome_print(MigOutcome *self, FILE *fp);

#endif
