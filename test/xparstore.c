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
    assert(ParStore_nFixed(ps) == 0);
    assert(ParStore_nFree(ps) == 0);

    double val, *ptr, *vec;

    val = 12.3;
    ParStore_addFixedPar(ps, val, "x");
    ptr = ParStore_findPtr(ps, "x");
    assert(*ptr == val);
    assert(ParStore_nFixed(ps) == 1);
    assert(ParStore_nFree(ps) == 0);
    assert(ParStore_getFixed(ps, 0) == val);

    val = 23.4;
    ParStore_addFreePar(ps, val, 10.0, 30.0, "y");
    ptr = ParStore_findPtr(ps, "y");
    assert(*ptr == val);
    assert(ParStore_nFixed(ps) == 1);
    assert(ParStore_nFree(ps) == 1);
    assert(ParStore_getFree(ps, 0) == val);
    assert(ParStore_loFree(ps, 0) == 10.0);
    assert(ParStore_hiFree(ps, 0) == 30.0);
    vec = ParStore_rawArray(ps);
    assert(vec[0] == val);

    val = 88.3;
    ParStore_addFixedPar(ps, val, "w");
    ptr = ParStore_findPtr(ps, "w");
    assert(*ptr == val);
    assert(ParStore_nFixed(ps) == 2);
    assert(ParStore_nFree(ps) == 1);
    assert(ParStore_getFixed(ps, 1) == val);

    val = -23.8;
    ParStore_addFreePar(ps, val, -100.0, 0.0, "z");
    ptr = ParStore_findPtr(ps, "z");
    assert(*ptr == val);
    assert(ParStore_nFixed(ps) == 2);
    assert(ParStore_nFree(ps) == 2);
    assert(ParStore_getFree(ps, 1) == val);
    assert(ParStore_loFree(ps, 1) == -100.0);
    assert(ParStore_hiFree(ps, 1) == 0.0);
    vec = ParStore_rawArray(ps);
    assert(vec[1] == val);

    ParStore_free(ps);

    unitTstResult("ParStore", "OK");

    return 0;
}
