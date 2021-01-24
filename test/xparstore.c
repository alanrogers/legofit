/**
 * @file xparstore.c
 * @author Alan R. Rogers
 * @brief Test parstore.c.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "parstore.h"
#include "param.h"
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
    double fixed0=99.0, fixed1=100.0;
    PtrQueue *fixedQ = PtrQueue_new();
    PtrQueue *freeQ = PtrQueue_new();
    PtrQueue *constrQ = PtrQueue_new();

    Param *par = Param_new("x", x, 0.0, 100.0, TWON|FREE, NULL);
    PtrQueue_push(freeQ, par);
    
    par = Param_new("y", y, 0.0, 100.0, TWON|FREE, NULL);
    PtrQueue_push(freeQ, par);

    par = Param_new("z", z, 0.0, 100.0, TIME|FREE, NULL);
    PtrQueue_push(freeQ, par);

    par = Param_new("fixed0", fixed0, fixed0, fixed0, TIME|FIXED, NULL);
    PtrQueue_push(fixedQ, par);

    par = Param_new("fixed1", fixed1, fixed1, fixed1, TWON|FIXED, NULL);
    PtrQueue_push(fixedQ, par);

    par = Param_new("c", 0.0, -DBL_MAX, DBL_MAX, TIME|CONSTRAINED,
                    "exp(x+log(y+z))");
    PtrQueue_push(constrQ, par);

    ParStore *ps = ParStore_new(fixedQ, freeQ, constrQ);

    int ndx = ParStore_getIndex(ps, "c");
    assert(ndx >= 0);
    assert(ndx < ParStore_nPar(ps));
    assert(exp(x+log(y+z)) == ParStore_getVal(ps, ndx));

    assert(6 == ParStore_nPar(ps));
    assert(2 == ParStore_nFixed(ps));
    assert(3 == ParStore_nFree(ps));
    assert(1 == ParStore_nConstrained(ps));

    const int n=3;
    double v[n];

    ParStore_getFreeParams(ps, n, v);
    assert(v[0] == x);
    assert(v[1] == y);
    assert(v[2] == z);

    v[0]=x=3.0;
    v[1]=y=2.0;
    v[2]=z=1.0;
    ParStore_setFreeParams(ps, n, v);

    ndx = ParStore_getIndex(ps, "x");
    assert(x == ParStore_getVal(ps, ndx));

    ndx = ParStore_getIndex(ps, "y");
    assert(y == ParStore_getVal(ps, ndx));

    ndx = ParStore_getIndex(ps, "z");
    assert(z == ParStore_getVal(ps, ndx));

    ndx = ParStore_getIndex(ps, "c");
    assert(exp(x+log(y+z)) == ParStore_getVal(ps, ndx));

    if(verbose)
        ParStore_print(ps, stdout);

    fflush(stdout);
    ParStore *ps2 = ParStore_dup(ps);

    assert(ParStore_equals(ps, ps2));

    assert(ParStore_nFree(ps2) == ParStore_nFree(ps));
    assert(ParStore_nFixed(ps2) == ParStore_nFixed(ps));

    ParStore_free(ps);
    ParStore_free(ps2);

    assert(0 == PtrQueue_size(freeQ));
    assert(0 == PtrQueue_size(fixedQ));
    assert(0 == PtrQueue_size(constrQ));

    PtrQueue_free(freeQ);
    PtrQueue_free(fixedQ);
    PtrQueue_free(constrQ);
    te_free_func_map();

    gsl_rng_free(rng);
    
    unitTstResult("ParStore", "OK");

    return 0;
}
