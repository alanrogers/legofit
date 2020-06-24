#include "migoutcome.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>

// not thread safe
static unsigned migration_event = 0;

static int MigOutcome_cmp(MigOutcome *a, unsigned event, uint64_t bits);
static MigOutcome *MigOutcome_new(MigOutcome *next,
                                  unsigned event,
                                  unsigned noutcomes,  
                                  unsigned outcome,
                                  double pr);

/// Increment external migration_event variable.
unsigned nextMigrationEvent(void) {
    unsigned event = migration_event;
    migration_event += 1;
    return event;
}

static int MigOutcome_cmp(MigOutcome *a, unsigned event, uint64_t bits) {
    if(a == NULL || a->event > event)
        return 1;
    if(a->event < event)
        return -1;
    if(a->bits > bits)
        return 1;
    if(a->bits < bits)
        return -1;
    return 0;
}

MigOutcome *MigOutcome_insert(MigOutcome *head,
                              unsigned event,
                              unsigned noutcomes,  
                              unsigned outcome,
                              double pr) {
    int cmp = MigOutcome_cmp(head, event,
                             ((uint64_t) 1) << outcome);
    if(cmp < 0) {
        head->next = MigOutcome_insert(head->next, event, noutcomes,
                                       outcome, pr);
        return head;
    }else if(cmp > 0) {
        return MigOutcome_new(head, event, noutcomes, outcome, pr);
    }else{
        fprintf(stderr,"%s:%d: attempt to insert duplicate MigOutcome event\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    return NULL; // NOTREACHED
}

static MigOutcome *MigOutcome_new(MigOutcome *next,
                                  unsigned event,
                                  unsigned noutcomes,  
                                  unsigned outcome,
                                  double pr) {
    MigOutcome *self = malloc(sizeof(MigOutcome));
    CHECKMEM(self);

    self->event = event;
    self->noutcomes = noutcomes;

    // set bit a position determined by outcome.
    static const unsigned nbits = 8*sizeof(self->bits);
    if(outcome >= nbits) {
        // The number of outcomes is 2^k, where k is the number of
        // lineages, and 2^6 = 64. So with type uint64_t, we can only
        // handle migration in segments with 6 or fewer lineages.
        // To handle 7, use type __uint128_t.
        fprintf(stderr,"%s:%d: can't set bit %d in a %d-bit integer.\n",
                __FILE__,__LINE__, outcome+1, nbits);
        fprintf(stderr,"%s:%d:"
                " This happens when the number of lineages is too large"
                " during\n     a migration event.\n",
                __FILE__,__LINE__);
        return NULL;
    }
    self->bits = 1;
    self->bits <<= outcome;
    self->pr = pr;
    self->next = next;
    return self;
}

MigOutcome *MigOutcome_dup(MigOutcome *old) {
    if(old == NULL)
        return NULL;
    MigOutcome *new = memdup(old, sizeof(MigOutcome));
    CHECKMEM(new);
    new->next = MigOutcome_dup(old->next);
    return new;
}

void MigOutcome_free(MigOutcome *self) {
    if(self == NULL)
        return;
    MigOutcome_free(self->next);
    free(self);
}

void MigOutcome_print(MigOutcome *self, FILE *fp) {
    if(self==NULL) {
        putc('\n', fp);
        return;
    }
    fprintf(fp,"%u.%llx:%g ",self->event, self->bits, self->pr);
    MigOutcome_print(self->next, fp);
}

int    IdSet_copyMigOutcome(IdSet *self, const IdSet *old);
