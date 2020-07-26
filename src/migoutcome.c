#include "migoutcome.h"
#include "misc.h"
#include "ptrqueue.h"
#include <stdio.h>
#include <stdlib.h>

// not thread safe
static unsigned migration_event = 0;

static int MigOutcome_cmp(MigOutcome *a, unsigned event);
static MigOutcome *MigOutcome_new(MigOutcome *next,
                                  unsigned event,
                                  unsigned noutcomes,  
                                  uint64_t bits,
                                  double pr);

/// Increment external migration_event variable.
unsigned nextMigrationEvent(void) {
    unsigned event = migration_event;
    migration_event += 1;
    return event;
}

static int MigOutcome_cmp(MigOutcome *a, unsigned event) {
    if(a == NULL || a->event > event)
        return 1;
    if(a->event < event)
        return -1;
    return 0;
}

MigOutcome *MigOutcome_insert(MigOutcome *head,
                              unsigned event,
                              unsigned noutcomes,  
                              uint64_t bits,
                              double pr) {
    int cmp = MigOutcome_cmp(head, event);
    if(cmp < 0) {
        head->next = MigOutcome_insert(head->next, event, noutcomes,
                                       bits, pr);
        return head;
    }else if(cmp > 0) {
        return MigOutcome_new(head, event, noutcomes, bits, pr);
    }else{
        fprintf(stderr,"%s:%s:%d: can't insert a new outcome for"
                " existing migration event %u\n",
                __FILE__,__func__,__LINE__, event);
        exit(EXIT_FAILURE);
    }

    return NULL; // NOTREACHED
}

static MigOutcome *MigOutcome_new(MigOutcome *next,
                                  unsigned event,
                                  unsigned noutcomes,  
                                  uint64_t bits,
                                  double pr) {
    MigOutcome *self = malloc(sizeof(MigOutcome));
    CHECKMEM(self);

    self->event = event;
    self->noutcomes = noutcomes;
    self->bits = bits;
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
    fprintf(fp,"%u:%llx:%g ",self->event, self->bits, self->pr);
    MigOutcome_print(self->next, fp);
}

/**
 * Return a newly-allocated MigOutcome list that includes all the
 * migration events in "left" and "right". If "left" and "right" are
 * mutually exclusive, return NULL.
 */
MigOutcome *MigOutcome_join(MigOutcome *left, MigOutcome *right) {

    MigOutcome *head = NULL;

    while(left || right) {
        if(right==NULL || (left->event < right->event)) {
            head = MigOutcome_insert(head, left->event,
                                     left->noutcomes, left->bits,
                                     left->pr);
            left = left->next;
        }else if(left==NULL || (left->event > right->event)) {
            head = MigOutcome_insert(head, right->event,
                                     right->noutcomes, right->bits,
                                     right->pr);
            right = right->next;
        }else{
            assert(left->event == right->event);
            if(left->outcome != right->outcome) {
                // left and right represent mutually exclusive events
                MigOutcomee_free(head);
                return NULL;
            }
            assert(left->pr == right->pr);
            head = MigOutcome_insert(head, right->event,
                                     right->noutcomes, right->bits,
                                     right->pr);
            left = left->next;
            right = right->next;
        }
    }
    return head;
}
