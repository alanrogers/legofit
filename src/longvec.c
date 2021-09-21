/**
@file longvec.c
@page LongVec
@brief An expandable vector of long ints.

@copyright Copyright (c) 2021, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
**/

#include "longvec.h"
#include "misc.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

struct LongVec {
    int size;
    long *vec;  // array of pointers
};

LongVec *LongVec_new(int size) {
    LongVec *new = malloc(sizeof(LongVec));
    CHECKMEM(new);

    new->size = size;
    new->vec = malloc(size * sizeof(new->vec[0]));
    CHECKMEM(new->vec);
    memset(new->vec, 0, size * sizeof(new->vec[0]) );
    return new;
}

int LongVec_size(const LongVec *self) {
    return self->size;
}

int LongVec_resize(LongVec *self, int newsize) {
    long *old = self->vec;
    int oldsize = self->size;
    self->size = newsize;
    self->vec = realloc(self->vec, newsize * sizeof(self->vec[0]));
    if(self->vec == NULL) {
        self->vec = old;
        return ENOMEM;
    }
    memset(self->vec + oldsize, 0, (newsize-oldsize)*sizeof(self->vec[0]));
    return 0;
}

void LongVec_set(LongVec *self, int ndx, long value) {
    assert(ndx < self->size);
    self->vec[ndx] = value;
}

void LongVec_plusEquals(LongVec *self, int ndx, long inc) {
    assert(ndx < self->size);
    self->vec[ndx] += inc;
}

void LongVec_minusEquals(LongVec *self, int ndx, long dec) {
    assert(ndx < self->size);
    self->vec[ndx] -= dec;
}

long LongVec_get(const LongVec *self, int ndx) {
    assert(ndx < self->size);
    return self->vec[ndx];
}

void LongVec_free(LongVec *self) {
    free(self->vec);
    free(self);
}



