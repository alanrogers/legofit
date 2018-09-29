/**
 * @file xparstore.c
 * @author Alan R. Rogers
 * @brief Test parstore.c.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "parstore.h"
#include "strparmap.h"
#include "addrparmap.h"
#include "ptrset.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <gsl/gsl_rng.h>
#include <time.h>

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

    time_t      currtime = time(NULL);
    unsigned    seed = currtime % UINT_MAX;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    CHECKMEM(rng);
    gsl_rng_set(rng, seed);

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

    double x=1.0, y=1.0, z=2.0;

    char buff[100];
    strcpy(buff, "1+1*x + 2*x*y");


    ParStore *ps = ParStore_new();
    assert(ParStore_nFixed(ps) == 0);
    assert(ParStore_nFree(ps) == 0);
    assert(ParStore_nConstrained(ps) == 0);
    ParStore_free(ps);

    double val, *ptr;
    unsigned ptype;

    PtrSet *seen = PtrSet_new();

    ps = ParStore_new();
    double fixed0=99.0, fixed1=100.0;
    ParStore_addFreePar(ps, x, 0.0, 100.0, "x", TWON);
    ParStore_addFreePar(ps, y, 0.0, 100.0, "y", TWON);
    ParStore_addFreePar(ps, z, 0.0, 100.0, "z", TIME);
    ParStore_addFixedPar(ps, fixed0, "fixed0", TIME);
    ParStore_addFixedPar(ps, fixed1, "fixed1", TWON);
    ParStore_addConstrainedPar(ps, "exp(x+log(y+z))", "c", TIME);
    ParStore_constrain(ps);
    ptr = ParStore_findPtr(ps, &ptype, "c");
    assert(ptype & CONSTRAINED);
    assert(ptype & TIME);
    assert(ParStore_isConstrained(ps, ptr));
    PtrSet_insert(seen, ptr);
    ParStore_chkDependencies(ps, ptr, seen);
    assert(2 == ParStore_nFixed(ps));
    assert(3 == ParStore_nFree(ps));
    assert(1 == ParStore_nConstrained(ps));
    assert(exp(x+log(y+z)) == *ptr);
    ParStore_addConstrainedPar(ps, "x+c", "d", ARBITRARY);
    ParStore_constrain(ps);
    ptr = ParStore_findPtr(ps, &ptype, "d");
    assert(ptype & CONSTRAINED);
    assert(ptype & ARBITRARY);
    assert(ParStore_isConstrained(ps, ptr));
    ParStore_chkDependencies(ps, ptr, seen);
    ParStore_free(ps);

    ps = ParStore_new();
    val = 12.3;
    ParStore_addFixedPar(ps, val, "x", TWON);
    ptr = ParStore_findPtr(ps, &ptype, "x");
    assert(*ptr == val);
    assert(!ParStore_isConstrained(ps, ptr));
	assert(ptype & FIXED);
	assert(ptype & TWON);
    assert(ParStore_nFixed(ps) == 1);
    assert(ParStore_nFree(ps) == 0);

    val = 23.4;
    ParStore_addFreePar(ps, val, 10.0, 30.0, "y", TIME);
    ptr = ParStore_findPtr(ps, &ptype, "y");
    assert(*ptr == val);
	assert(ptype & FREE);
	assert(ptype & TIME);
    assert(ParStore_nFixed(ps) == 1);
    assert(ParStore_nFree(ps) == 1);

    val = 88.3;
    ParStore_addFixedPar(ps, val, "w", TWON);
    ptr = ParStore_findPtr(ps, &ptype, "w");
    assert(*ptr == val);
    assert(!ParStore_isConstrained(ps, ptr));
	assert(ptype & FIXED);
	assert(ptype & TWON);
    assert(ParStore_nFixed(ps) == 2);
    assert(ParStore_nFree(ps) == 1);

    val = -23.8;
    ParStore_addFreePar(ps, val, -100.0, 0.0, "z", ARBITRARY);
    ptr = ParStore_findPtr(ps, &ptype, "z");
    assert(*ptr == val);
    assert(!ParStore_isConstrained(ps, ptr));
	assert(ptype & FREE);
	assert(ptype & ARBITRARY);
    assert(ParStore_nFixed(ps) == 2);
    assert(ParStore_nFree(ps) == 2);

    val = 0.8;
    ParStore_addFreePar(ps, val, 0.0, 1.0, "a", MIXFRAC);
    ptr = ParStore_findPtr(ps, &ptype, "a");
    assert(*ptr == val);
    assert(!ParStore_isConstrained(ps, ptr));
	assert(ptype & FREE);
	assert(ptype & MIXFRAC);
    assert(ParStore_nFixed(ps) == 2);
    assert(ParStore_nFree(ps) == 3);

    strcpy(buff, "1 + 2*y + 1*a");
    assert(0 == ParStore_nConstrained(ps));
    ParStore_addConstrainedPar(ps, buff, "cnstr", ARBITRARY);
    assert(1 == ParStore_nConstrained(ps));
    ParStore_constrain(ps);

    if(verbose)
        ParStore_print(ps, stdout);

    fflush(stdout);
    ParStore *ps2 = ParStore_dup(ps);

    assert(ParStore_equals(ps, ps2));

    ParStore_randomize(ps2, rng);
    ParStore_sanityCheck(ps2, __FILE__, __LINE__);

    assert( !ParStore_equals(ps, ps2) );

    assert(ParStore_nFree(ps2) == ParStore_nFree(ps));
    assert(ParStore_nFixed(ps2) == ParStore_nFixed(ps));

    ParStore_free(ps);
    ParStore_free(ps2);
    unitTstResult("ParStore", "OK");

    return 0;
}
