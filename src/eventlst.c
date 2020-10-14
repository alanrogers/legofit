#include "eventlst.h"
#include "misc.h"
#include "ptrqueue.h"
#include <stdio.h>
#include <stdlib.h>

static int EventLst_cmp_shallow(EventLst *a, unsigned event);
static EventLst *EventLst_new(EventLst *next,
                                  unsigned event,
                                  unsigned outcomes,  
                                  long double pr);

/// Return 0-based index of next event. Not thread safe.
unsigned nextEvent(void) {
    static unsigned event = 0;
    return event++;
}

static int EventLst_cmp_shallow(EventLst *a, unsigned event) {
    if(a == NULL || a->event < event)
        return 1;
    if(a->event > event)
        return -1;
    return 0;
}

EventLst *EventLst_insert(EventLst *head,
                    unsigned event,
                    unsigned outcome,
                    long double pr) {
    int cmp = EventLst_cmp_shallow(head, event);
    if(cmp < 0) {
        head->next = EventLst_insert(head->next, event, outcome,
                                  pr);
        return head;
    }else if(cmp > 0) {
        return EventLst_new(head, event, outcome, pr);
    }else{
        fprintf(stderr,"%s:%s:%d: can't insert a new outcome for"
                " existing event %u\n",
                __FILE__,__func__,__LINE__, event);
        exit(EXIT_FAILURE);
    }

    return NULL; // NOTREACHED
}

static EventLst *EventLst_new(EventLst *next,
                                  unsigned event,
                                  unsigned outcome,  
                                  long double pr) {
    EventLst *self = malloc(sizeof(EventLst));
    CHECKMEM(self);

    self->event = event;
    self->outcome = outcome;
    self->pr = pr;
    self->next = next;
    return self;
}

EventLst *EventLst_dup(EventLst *old) {
    if(old == NULL)
        return NULL;
    EventLst *new = memdup(old, sizeof(EventLst));
    CHECKMEM(new);
    new->next = EventLst_dup(old->next);
    return new;
}

void EventLst_free(EventLst *self) {
    if(self == NULL)
        return;
    EventLst_free(self->next);
    free(self);
}

void EventLst_print(EventLst *self, FILE *fp) {
    for(EventLst *m=self; m; m = m->next) {
        if(m != self)
            putc(' ', fp);
        fprintf(fp,"%u_%u_%Lg", m->event, m->outcome, m->pr);
    }
}

/**
 * Return a newly-allocated EventLst list that includes all the events in
 * "left" and "right". If "left" and "right" are mutually exclusive,
 * return NULL.
 */
EventLst *EventLst_join(EventLst *left, EventLst *right,
                            int *mutually_exclusive) {

    EventLst *head = NULL;

    while(left && right) {
        if(left->event < right->event) {
            head = EventLst_insert(head, left->event,
                                     left->outcome, left->pr);
            left = left->next;
        }else if(left->event > right->event) {
            head = EventLst_insert(head, right->event,
                                     right->outcome, right->pr);
            right = right->next;
        }else{
            assert(left->event == right->event);
            if(left->outcome != right->outcome) {
                // left and right represent mutually exclusive events
                EventLst_free(head);
                *mutually_exclusive = 1;
                return NULL;
            }
            assert(left->pr == right->pr);

            // This is the same event, not two independent events,
            // so the probabilities of left and right do not
            // multiply.
            head = EventLst_insert(head, right->event,
                                     right->outcome, right->pr);
            left = left->next;
            right = right->next;
        }
    }
    while(left) {
       head = EventLst_insert(head, left->event, left->outcome,
                                left->pr);
       left = left->next;
    }
    while(right) {
        head = EventLst_insert(head, right->event, right->outcome,
                                 right->pr);
        right = right->next;
    }
    *mutually_exclusive = 0;
    return head;
}

/// Probability of current EventLst list is the product of
/// the probabilties of the elements of the list.
long double EventLst_prob(EventLst *head) {
    EventLst *mo;
    long double pr = 1.0L;

    for(mo=head; mo; mo = mo->next)
        pr *= mo->pr;

    return pr;
}

int EventLst_cmp(const EventLst *left, const EventLst *right) {
    if(left==NULL && right==NULL)
        return 0;
    if(right==NULL) // left!=NULL
        return 1;
    if(left==NULL)  // right!=NULL
        return -1;
    if(left->event < right->event)
        return 1;
    if (left->event > right->event)
        return -1;
    if (left->outcome < right->outcome)
        return 1;
    if (left->outcome > right->outcome)
        return -1;
    return EventLst_cmp(left->next, right->next);
}

uint32_t EventLst_hash(const EventLst *self)
    __attribute__((no_sanitize("integer"))) {
    if(self==NULL)
        return 0;
    uint32_t hash = 17;
    hash = hash * 37 + uint32Hash(self->event);
    hash = hash * 37 + uint32Hash(self->outcome);
    uint32_t h = EventLst_hash(self->next);
    if(h)
        hash = hash * 37 + h;
    return hash;
}
