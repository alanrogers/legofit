/**
@file uintqueue.c
@brief A FIFO queue of unsigned integers.
*/

#include "uintqueue.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>

struct UINTqueue {
    struct UINTqueue *next;
    unsigned value;
};

// Push a value onto the tail of the queue. Return pointer to new
// head. Example:
// 
// UINTqueue *queue=NULL;
// queue = UINTqueue_push(queue, 1u);
// queue = UINTqueue_push(queue, 2u);
UINTqueue *UINTqueue_push(UINTqueue *self, unsigned value) {
    if(self != NULL) {
        self->next = UINTqueue_push(self->next, value);
        return self;
    }
    UINTqueue *new = malloc(sizeof(UINTqueue));
    CHECKMEM(new);
    new->value = value;
    new->next = NULL;
    return new;
}

// Pop a value off the head of the queue. Return pointer to new
// head. Example:
// 
// UINTqueue *queue=NULL;
// queue = UINTqueue_push(queue, 1u);
// queue = UINTqueue_push(queue, 2u);
//
// unsigned x;
// queue = UINTqueue_pop(queue, &x);  // x=1
// queue = UINTqueue_pop(queue, &x);  // x=2
UINTqueue *UINTqueue_pop(UINTqueue *self, unsigned *value) {
    if(self==NULL)
        return NULL;
    *value = self->value;
    UINTqueue *next = self->next;
    free(self);
    return next;
}

int UINTqueue_length(UINTqueue *self) {
    if(self==NULL)
        return 0;
    return 1 + UINTqueue_length(self->next);
}

UINTqueue *UINTqueue_free(UINTqueue *self) {
    if(self) {
        self->next = UINTqueue_free(self->next);
        free(self);
    }
    return NULL;
}


