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
#include <string.h>

#ifdef NDEBUG
#  error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {
	int verbose=0;

	switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xparstore [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xparstore [-v]\n");
    }

    Bounds bnd0 = {
        .lo_twoN = 0.0,
        .hi_twoN = 1e12,
        .lo_t = 0,
        .hi_t = 1e10
    };
    Bounds_sanityCheck(&bnd0, __FILE__, __LINE__);
    Bounds bnd1 = bnd0;
    Bounds_sanityCheck(&bnd1, __FILE__, __LINE__);
    assert(Bounds_equals(&bnd0, &bnd1));
    unitTstResult("Bounds", "OK");

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

    if(verbose)
        ParStore_print(ps, stdout);

    ParStore *ps2 = ParStore_dup(ps);
    size_t offset = ((size_t) ps2) - ((size_t) ps);
    int    i;

    assert(ParStore_equals(ps, ps2));

    assert(ParStore_nFree(ps2) == ParStore_nFree(ps));
    assert(ParStore_nFixed(ps2) == ParStore_nFixed(ps));

    for(i=0; i < ParStore_nFixed(ps); ++i) {
        assert(ParStore_getFixed(ps2, i) == ParStore_getFixed(ps, i));
        const char *name1 = ParStore_getNameFixed(ps, i);
        const char *name2 = ParStore_getNameFixed(ps2, i);
        assert(0 == strcmp(name1, name2));
        size_t position1 = (size_t) ParStore_findPtr(ps, name1);
        size_t position2 = (size_t) ParStore_findPtr(ps2, name1);
        assert(offset == position2 - position1);
    }

    for(i=0; i < ParStore_nFree(ps); ++i) {
        assert(ParStore_getFree(ps2, i) == ParStore_getFree(ps, i));
        assert(ParStore_loFree(ps2, i) == ParStore_loFree(ps, i));
        assert(ParStore_hiFree(ps2, i) == ParStore_hiFree(ps, i));
        const char *name1 = ParStore_getNameFree(ps, i);
        const char *name2 = ParStore_getNameFree(ps2, i);
        assert(0 == strcmp(name1, name2));
        size_t position1 = (size_t) ParStore_findPtr(ps, name1);
        size_t position2 = (size_t) ParStore_findPtr(ps2, name1);
        assert(offset == position2 - position1);
    }
    if(verbose)
        ParStore_print(ps2, stdout);

    ParStore_free(ps);
    ParStore_free(ps2);

    unitTstResult("ParStore", "OK");

    return 0;
}
