/**
 * @file xparstore.c
 * @author Alan R. Rogers
 * @brief Test parstore.c.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "parstore.h"
#include <stdio.h>
#include <assert.h>

#ifdef NDEBUG
#  error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(void) {
    ParStore *ps = ParStore_new();
    assert(ParStore_nPar(ps) == 0);

    double *x, *y;

    x = ParStore_addPar(ps, 12.3, 2.0, 15.0);
    assert(ParStore_nPar(ps) == 1);
    assert(ParStore_get(ps, 0) == *x);
    assert(ParStore_get(ps, 0) == 12.3);
    assert(ParStore_loBnd(ps, 0) == 2.0);
    assert(ParStore_hiBnd(ps, 0) == 15.0);

    y = ParStore_addPar(ps, -0.23, -1.0, 0.0);
    assert(ParStore_nPar(ps) == 2);
    assert(ParStore_get(ps, 1) == *y);
    assert(ParStore_get(ps, 1) == -0.23);
    assert(ParStore_loBnd(ps, 1) == -1.0);
    assert(ParStore_hiBnd(ps, 1) == 0.0);

    x = ParStore_getPtr(ps);
    assert(x[0] == 12.3);
    assert(x[1] == -0.23);

    ParStore_free(ps);

    unitTstResult("ParStore", "OK");

    return 0;
}
