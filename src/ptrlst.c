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

struct PtrLst {
    long unsigned count;
    PtrLstEl *head;
};

/// A linked list of pointers
struct PtrLstEl {
    struct PtrLstEl  *next;          ///< next item in list
    void *ptr;
};

static PtrLstEl *PtrLstEl_push(PtrLstEl *next, void * ptr);
static PtrLstEl *PtrLstEl_pop(PtrLstEl *head, void **ptr);
static void PtrLstEl_free(PtrLstEl * self);

/// Push pointer onto list.
static PtrLstEl *PtrLstEl_push(PtrLstEl *next, void * ptr) {
    PtrLstEl *new = malloc(sizeof(*new));
    CHECKMEM(new);

    new->next = next;
    new->ptr = ptr;
    return new;
}

/// Returns NULL if list is empty, so this list cannot be used
/// to store NULL pointers.
static PtrLstEl *PtrLstEl_pop(PtrLstEl *head, void **ptr) {
    if(head == NULL) {
        *ptr = NULL;
        return NULL;
    }
    *ptr = head->ptr;
    PtrLstEl *next = head->next;
    free(head);
    return next;
}

/// Destroy a linked list of PtrLstEl objects.
static void PtrLstEl_free(PtrLstEl * e) {
    if(e == NULL)
        return;
    PtrLstEl_free(e->next);
    free(e);
}

/// Allocate a structure to hold a linked list.
PtrLst *PtrLst_new(void) {
    PtrLst *self = malloc(sizeof(PtrLst));
    CHECKMEM(self);
    self->head = NULL;
    self->count = 0;
    return self;
}

void PtrLst_free(PtrLst *self) {
    PtrLstEl_free(self->head);
    free(self);
}

int PtrLst_push(PtrLst *self, void *ptr) {
    self->head = PtrLstEl_push(self->head, ptr);
    self->count += 1;
    return 0;
}

void *PtrLst_pop(PtrLst *self) {
    void *ptr;
    self->head = PtrLstEl_pop(self->head, &ptr);
    if(self->count > 0)
        self->count -= 1;
    return ptr;
}

long unsigned PtrLst_length(PtrLst *self) {
    return self->count;
}

void    PtrLst_rewind(PtrLst *self);
void   *PtrLst_next(PtrLst *self);
