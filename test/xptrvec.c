/**
 * @file xptrvec.c
 * @author Alan R. Rogers
 * @brief Test ptrvec.c.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include "ptrvec.h"
#include "misc.h"
#include <stdio.h>
#include <assert.h>
 
int main(void) {

    int a=1, b=2, c=3, d=4, status;
    PtrVec *v = PtrVec_new(0);
    assert(v->buffsize == 0);
    assert(v->used == 0);

    status = PtrVec_push(v, &a);
    assert(status == 0);
    assert(v->buffsize == 2);
    assert(v->used == 1);
    assert(PtrVec_get(v, 0) == &a);

    status = PtrVec_push(v, &b);
    assert(status == 0);
    assert(v->buffsize == 2);
    assert(v->used == 2);
    assert(PtrVec_get(v, 1) == &b);
    
    status = PtrVec_push(v, &c);
    assert(status == 0);
    assert(v->buffsize == 4);
    assert(v->used == 3);
    assert(PtrVec_get(v, 2) == &c);
    
    status = PtrVec_push(v, &d);
    assert(status == 0);
    assert(v->buffsize == 4);
    assert(v->used == 4);
    assert(PtrVec_get(v, 3) == &d);

    int *p = PtrVec_pop(v);
    assert(p == &d);
    assert(v->buffsize == 4);
    assert(v->used == 3);

    p = PtrVec_pop(v);
    assert(p == &c);
    assert(v->buffsize == 4);
    assert(v->used == 2);

    p = PtrVec_pop(v);
    assert(p == &b);
    assert(v->buffsize == 4);
    assert(v->used == 1);

    p = PtrVec_pop(v);
    assert(p == &a);
    assert(v->buffsize == 4);
    assert(v->used == 0);
    
    p = PtrVec_pop(v);
    assert(p == NULL);

    PtrVec_free(v);

    unitTstResult("PtrVec", "OK");
    return 0;
}
