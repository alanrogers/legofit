/**
 * @file ptrqueue.c
 * @author Alan R. Rogers
 * @brief A FIFO queue of pointers.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "ptrqueue.h"
#include "misc.h"
#include <stdlib.h>

typedef struct El El;

// One link in the chain
struct El {
    void *ptr; // not locally owned
    struct El *next;
};

// Head and tail of the linked list.
struct PtrQueue {
    El *head, *tail;
};

static El *El_push(El *tail, void *ptr);
static El *El_pop(El *head, void **ptr);
static void El_free(El *self);

static El *El_push(El *tail, void *ptr) {
    El *self = malloc(sizeof(El));
    CHECKMEM(self);
    self->ptr = ptr;
    self->next = NULL;
    if(tail)
        tail->next = self;
    return self;
}

static El *El_pop(El *head, void **ptr) {
    if(head==NULL) {
        *ptr = NULL;
        return NULL;
    }
    *ptr = head->ptr;
    El *next = head->next;
    free(head);
    return next;
}

static void El_free(El *self) {
    for(El *e = self; e; e=e->next)
        free(e);
}

PtrQueue *PtrQueue_new(void) {
    PtrQueue *self = malloc(sizeof(PtrQueue));
    CHECKMEM(self);
    self->head = self->tail = NULL;
    return self;
}

void PtrQueue_free(PtrQueue *self) {
    El_free(self->head);
    free(self);
}

void PtrQueue_push(PtrQueue *self, void *ptr) {
    self->tail = El_push(self->tail, ptr);
    if(self->head == NULL)
        self->head = self->tail;
}

void *PtrQueue_pop(PtrQueue *self) {
    void *ptr;
    self->head = El_pop(self->head, &ptr);
    if(self->head == NULL)
        self->tail = NULL;
    return ptr;
}

unsigned PtrQueue_size(PtrQueue *self) {
    unsigned n=0;
    for(El *e = self->head; e; e=e->next)
        n += 1;
    return n;
}
