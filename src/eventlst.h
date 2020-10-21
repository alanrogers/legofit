#ifndef ARR_EVENTLST
#define ARR_EVENTLST

#include "typedefs.h"
#include <stdio.h>
#include <stdint.h>

/**
This structure handles bookkeeping associated with events of two
types: (1) partitions of the set of descendants into a set of
ancestors, and (2) migration events. In either case, a number of
mutually exclusive outcomes are possible, and only ancestors of the
same outcome can be allowed to coalesce farther back in time.

Let k represent the number of ancestors at the ancient end of a
segment. Each of these is an immigrant with some positive
probability. There are 2^k mutually-exclusive ways to choose
immigrants and natives, and each of these outcomes divides the
ancestors into two sets, one of which may be empty.

In the EventLst structure, "event" is an index that indicates which
episode of migration is being described. It is 0 for the first
migration event, 1 for the next, and so on.

The "outcome" field is 0 for the first outcome, 1 for the second, and
so on.  It is forbidden to coalesce lineages that belong to different
outcomes of the same event.

A single set of lineages may have experienced multiple events, so each
lineage has a sorted linked list of EventLst objects, which represent all
the migration events it has experienced.
 **/
struct EventLst {
    unsigned event;        // which episode of migration?
    unsigned outcome;     // outcome of this event
    long double pr;       // probability of this outcome
    struct EventLst *next;
};

unsigned    nextEvent(void);
EventLst   *EventLst_insert(EventLst *head,
                            unsigned event,
                            unsigned outcome,  
                            long double pr);
EventLst   *EventLst_dup(EventLst *old);
void        EventLst_free(EventLst *self);
void        EventLst_print(EventLst *self, FILE *fp);
EventLst   *EventLst_join(EventLst *left, EventLst *right,
                          int *mutually_exclusive);
long double EventLst_prob(EventLst *head);
int         EventLst_cmp(const EventLst *left, const EventLst *right);
uint32_t    EventLst_hash(const EventLst *self)
    __attribute__((no_sanitize("integer")));
void        EventLst_sanityCheck(const EventLst *self,
                                 const char *file, int lineno);

#endif
