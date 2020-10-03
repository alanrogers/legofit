#include "migoutcome.h"
#include "misc.h"
#include "ptrqueue.h"
#include <stdio.h>
#include <stdlib.h>

// not thread safe
static unsigned migration_event = 0;

static int MigOutcome_cmp_shallow(MigOutcome *a, unsigned event);
static MigOutcome *MigOutcome_new(MigOutcome *next,
                                  unsigned event,
                                  unsigned outcomes,  
                                  long double pr);

/// Increment external migration_event variable.
unsigned nextMigrationEvent(void) {
    return migration_event++;
}

static int MigOutcome_cmp_shallow(MigOutcome *a, unsigned event) {
    if(a == NULL || a->event > event)
        return 1;
    if(a->event < event)
        return -1;
    return 0;
}

MigOutcome *MigOutcome_insert(MigOutcome *head,
                              unsigned event,
                              unsigned outcome,
                              long double pr) {
    int cmp = MigOutcome_cmp_shallow(head, event);
    if(cmp < 0) {
        head->next = MigOutcome_insert(head->next, event, outcome,
                                       pr);
        return head;
    }else if(cmp > 0) {
        return MigOutcome_new(head, event, outcome, pr);
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
                                  unsigned outcome,  
                                  long double pr) {
    MigOutcome *self = malloc(sizeof(MigOutcome));
    CHECKMEM(self);

    self->event = event;
    self->outcome = outcome;
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
    for(MigOutcome *m=self; m; m = m->next) {
        if(m != self)
            putc(' ', fp);
        fprintf(fp,"%u_%u_%Lg", m->event, m->outcome, m->pr);
    }
}

/**
 * Return a newly-allocated MigOutcome list that includes all the
 * migration events in "left" and "right". If "left" and "right" are
 * mutually exclusive, return NULL.
 */
MigOutcome *MigOutcome_join(MigOutcome *left, MigOutcome *right,
                            int *mutually_exclusive) {

    MigOutcome *head = NULL;

    while(left && right) {
        if(left->event < right->event) {
            head = MigOutcome_insert(head, left->event,
                                     left->outcome, left->pr);
            left = left->next;
        }else if(left->event > right->event) {
            head = MigOutcome_insert(head, right->event,
                                     right->outcome, right->pr);
            right = right->next;
        }else{
            assert(left->event == right->event);
            if(left->outcome != right->outcome) {
                // left and right represent mutually exclusive events
                MigOutcome_free(head);
                *mutually_exclusive = 1;
                return NULL;
            }
            assert(left->pr == right->pr);

            // This is the same event, not two independent events,
            // so the probabilities of left and right do not
            // multiply.
            head = MigOutcome_insert(head, right->event,
                                     right->outcome, right->pr);
            left = left->next;
            right = right->next;
        }
    }
    while(left) {
       head = MigOutcome_insert(head, left->event, left->outcome,
                                left->pr);
       left = left->next;
    }
    while(right) {
        head = MigOutcome_insert(head, right->event, right->outcome,
                                 right->pr);
        right = right->next;
    }
    *mutually_exclusive = 0;
    return head;
}

/// Probability of current MigOutcome list is the product of
/// the probabilties of the elements of the list.
long double MigOutcome_prob(MigOutcome *head) {
    MigOutcome *mo;
    long double pr = 1.0L;

    for(mo=head; mo; mo = mo->next)
        pr *= mo->pr;

    return pr;
}

int MigOutcome_cmp(const MigOutcome *left, const MigOutcome *right) {
    if(left==NULL && right==NULL)
        return 0;
    if(right==NULL) // left!=NULL
        return 1;
    if(left==NULL)  // right!=NULL
        return -1;
    if(left->event > right->event)
        return 1;
    if (left->event < right->event)
        return -1;
    if (left->outcome > right->outcome)
        return 1;
    if (left->outcome < right->outcome)
        return -1;
    return MigOutcome_cmp(left->next, right->next);
}

uint32_t MigOutcome_hash(const MigOutcome *self)
    __attribute__((no_sanitize("integer"))) {
    if(self==NULL)
        return 0;
    uint32_t hash = 17;
    hash = hash * 37 + uint32Hash(self->event);
    hash = hash * 37 + uint32Hash(self->outcome);
    uint32_t h = MigOutcome_hash(self->next);
    if(h)
        hash = hash * 37 + h;
    return hash;
}
