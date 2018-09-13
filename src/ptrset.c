/**
 * @file ptrset.c
 * @author Alan R. Rogers
 * @brief A set of pointers to doubles.
 *
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "ptrset.h"
#include <stdlib.h>
#include <assert.h>

typedef struct PtrSetNode {
    const double *ptr;
    struct PtrSetNode *next;
} PtrSetNode;

struct PtrSet {
    PtrSetNode *root;
};

static PtrSetNode *PtrSetNode_new(const double *ptr, PtrSetNode *next);
static PtrSetNode *PtrSetNode_insert(PtrSetNode * self, const double *ptr);
static int PtrSetNode_exists(PtrSetNode *self, const double *ptr);
static void PtrSetNode_free(PtrSetNode *self);

static PtrSetNode *PtrSetNode_new(const double *ptr, PtrSetNode *next) {
    PtrSetNode *self = malloc(sizeof(PtrSetNode));
    if(self == NULL) {
        fprintf(stderr,"%s:%d: bad malloc\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    self->ptr = ptr;
    self->next = next;
    return self;
}

static void PtrSetNode_free(PtrSetNode *self) {
    if(self == NULL)
        return;
    PtrSetNode_free(self->next);
    free(self);
}

/// Insert a new pointer. Do nothing if it already exists.
static PtrSetNode *PtrSetNode_insert(PtrSetNode * self, const double *ptr) {
    if(self == NULL)
        return PtrSetNode_new(ptr, NULL);
    if(ptr == self->ptr)
		return self;
    if(ptr > self->ptr) {
        self->next = PtrSetNode_insert(self->next, ptr);
        return self;
    }
    return PtrSetNode_new(ptr, self);
}

/// Return 1 if ptr is present in set; 0 otherwise.
static int PtrSetNode_exists(PtrSetNode *self, const double *ptr) {
    if(self == NULL)
        return 0;
    if(ptr == self->ptr)
        return 1;
    if(ptr > self->ptr)
        return PtrSetNode_exists(self->next, ptr);
    assert(ptr < self->ptr);
    return 0;
}

PtrSet *PtrSet_new(void) {
    PtrSet *self = malloc(sizeof(PtrSet));
    if(self == NULL) {
        fprintf(stderr,"%s:%d: bad malloc\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    self->root = NULL;
    return self;
}

void PtrSet_free(PtrSet *self) {
    PtrSetNode_free(self->root);
    free(self);
}

void PtrSet_insert(PtrSet *self, const double *ptr) {
    self->root = PtrSetNode_insert(self->root, ptr);
}

int PtrSet_exists(PtrSet *self, const double *ptr) {
    return PtrSetNode_exists(self->root, ptr);
}

void PtrSet_print(PtrSet *self, FILE *fp) {
    for(PtrSetNode *i = self->root; i; i = i->next)
        fprintf(fp, " %p", i->ptr);
    putc('\n', fp);
}

