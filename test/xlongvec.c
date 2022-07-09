/**
 * @file xlongvec.c
 * @author Alan R. Rogers
 * @brief Test longvec.c.
 * @copyright Copyright (c) 2021, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include "longvec.h"
#include "misc.h"
#include <stdio.h>
#include <assert.h>
 
int main(void) {

    int size = 2;

    LongVec *v = LongVec_new(size);

    assert(size == LongVec_size(v));
    assert(0 == LongVec_get(v, 0));
    assert(0 == LongVec_get(v, 1));
    LongVec_set(v, 0, 3);
    assert(3 == LongVec_get(v, 0));
    assert(0 == LongVec_get(v, 1));
    LongVec_plusEquals(v, 1, 1);
    assert(3 == LongVec_get(v, 0));
    assert(1 == LongVec_get(v, 1));
    LongVec_minusEquals(v, 0, 1);
    assert(2 == LongVec_get(v, 0));
    assert(1 == LongVec_get(v, 1));

    size += 1;
    LongVec_resize(v, size);
    assert(size == LongVec_size(v));
    assert(2 == LongVec_get(v, 0));
    assert(1 == LongVec_get(v, 1));
    assert(0 == LongVec_get(v, 2));

    LongVec_free(v);
    v = NULL;

    unitTstResult("LongVec", "OK");
    return 0;
}
