/**
@file ptrvec.c
@page PtrVec
@brief An expandable array (or stack) of pointers

@copyright Copyright (c) 2020, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
**/

#include "ptrvec.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>

/// Allocate an array with room for n pointers.
PtrVec *PtrVec_new(unsigned n) {
    PtrVec *self = malloc(sizeof(PtrVec));
    CHECKMEM(self);

    self->buffsize = n;
    self->used = 0;
    self->buff = malloc(n * sizeof(void *));
    CHECKMEM(self->buff);
    return self;
}

void PtrVec_free(PtrVec *self) {
    if(self->buff)
        free(self->buff);
    free(self);
}

/**
 * PtrVec doesn't own the things it stores pointers to. But if you
 * want free all the stored pointers anyway, this will do it.
 */
void PtrVec_freeHoldings(PtrVec *self) {
    for(unsigned i=0; i < self->used; ++i)
        free(self->buff[i]);
    self->used = 0;
}

int PtrVec_push(PtrVec *self, void *val) {
    if(self->used == self->buffsize) {
        // buffer overflow: reallocate
        unsigned n = self->buffsize;
        self->buffsize = (n==0? 2 : 2*n);
        void *old = self->buff;
        self->buff = realloc(self->buff, self->buffsize);
        if(self->buff == NULL) {
            if(old)
                free(old);
            return ENOMEM;
        }
    }
    self->buff[self->used++] = val;
    return 0;
}

void *PtrVec_pop(PtrVec *self) {
    if(self->used == 0)
        return NULL;
    self->used -= 1;
    return self->buff[self->used];
}

