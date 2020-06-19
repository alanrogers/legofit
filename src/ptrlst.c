/**
@file ptrlst.c
@page PtrLst
@brief An expandable array (or stack) of pointers

@copyright Copyright (c) 2020, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
**/

#include "ptrlst.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>

typedef struct El El;

struct PtrLst {
    El *head;
};


/// A linked list of pointers
struct El {
    struct El  *next;          ///< next item in list
    void *ptr;
};

static El *El_push(El *next, void * ptr);
static El *El_pop(El *head, void **ptr);
static void El_free(El * self);
static long unsigned El_length(El * self);

/// Push pointer onto list.
static El *El_push(El *next, void * ptr) {
    El *new = malloc(sizeof(*new));
    CHECKMEM(new);

    new->next = next;
    new->ptr = ptr;
    return new;
}

/// Returns NULL if list is empty, so this list cannot be used
/// to store NULL pointers.
static El *El_pop(El *head, void **ptr) {
    if(head == NULL) {
        *ptr = NULL;
        return NULL;
    }
    *ptr = head->ptr;
    El *next = head->next;
    free(head);
    return next;
}

/// Destroy a linked list of El objects.
static void El_free(El * e) {
    if(e == NULL)
        return;
    El_free(e->next);
    free(e);
}

static long unsigned El_length(El * self) {
    long unsigned n=0;
    for(El *el = self; el; el = el->next)
        n += 1;
    return n;
}

/// Allocate a structure to hold a linked list.
PtrLst *PtrLst_new(void) {
    PtrLst *self = malloc(sizeof(PtrLst));
    CHECKMEM(self);
    self->head = NULL;
    return self;
}

void PtrLst_free(PtrLst *self) {
    El_free(self->head);
    free(self);
}

int PtrLst_push(PtrLst *self, void *ptr) {
    self->head = El_push(self->head, ptr);
    return 0;
}

void *PtrLst_pop(PtrLst *self) {
    void *ptr;
    self->head = El_pop(self->head, &ptr);
    return ptr;
}

long unsigned PtrLst_length(PtrLst *self) {
    return El_length(self->head);
}
