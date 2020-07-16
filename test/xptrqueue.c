/**
 * @file xptrqueue.c
 * @author Alan R. Rogers
 * @brief Test ptrqueue.c.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "ptrqueue.h"
#include "misc.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(void) {

    int x=1, y=2, z=3;

    PtrQueue *pq = PtrQueue_new();

    assert(NULL == PtrQueue_pop(pq));
    assert(0 == PtrQueue_size(pq));

    PtrQueue_push(pq, &x);
    assert(1 == PtrQueue_size(pq));
    PtrQueue_push(pq, &y);
    assert(2 == PtrQueue_size(pq));
    PtrQueue_push(pq, &z);
    assert(3 == PtrQueue_size(pq));
    assert(&x == PtrQueue_pop(pq));
    assert(2 == PtrQueue_size(pq));
    assert(&y == PtrQueue_pop(pq));
    assert(1 == PtrQueue_size(pq));
    assert(&z == PtrQueue_pop(pq));
    assert(0 == PtrQueue_size(pq));
    assert(NULL == PtrQueue_pop(pq));

    PtrQueue_free(pq);

    unitTstResult("PtrQueue", "OK");

    return 0;
}
