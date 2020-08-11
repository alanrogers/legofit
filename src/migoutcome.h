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

In the MigOutcome structure, "event" is an index that indicates which
episode of migration is being described. It is 0 for the first
migration event, 1 for the next, and so on.

"Outcome" is 0 for the first outcome, 1 for the second, and so on.
It is forbidden to coalesce lineages that belong to different outcomes
of the same event. 

A single set of lineages may have experienced multiple migration
events, so each lineage has a sorted linked list of MigOutcome
objects, which represent all the migration events it has experienced.
 **/
struct MigOutcome {
    unsigned event;        // which episode of migration?
    unsigned outcome;     // outcome of this event
    double pr;             // probability of this outcome
    struct MigOutcome *next;
};

unsigned    nextMigrationEvent(void);
MigOutcome *MigOutcome_insert(MigOutcome *head,
                              unsigned event,
                              unsigned outcome,  
                              double pr);
MigOutcome *MigOutcome_dup(MigOutcome *old);
void        MigOutcome_free(MigOutcome *self);
void        MigOutcome_print(MigOutcome *self, FILE *fp);
MigOutcome *MigOutcome_join(MigOutcome *left, MigOutcome *right,
                            int *mutually_exclusive);

#endif
