/**
 * @file xparstore.c
 * @author Alan R. Rogers
 * @brief Test parstore.c.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "parstore.h"
#include "parkeyval.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <limits.h>

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
    unsigned    baseSeed = currtime % UINT_MAX;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    CHECKMEM(rng);
    gsl_rng_set(rng, baseSeed);

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
    ParKeyVal *pkv = NULL;
    pkv = ParKeyVal_add(pkv, "x", &x, Free);
    pkv = ParKeyVal_add(pkv, "y", &y, Free);
    pkv = ParKeyVal_add(pkv, "z", &z, Free);

    char buff[100];
    strcpy(buff, "1+1*x + 2*x*y");

    double xx=1.0, yy=1.0, zz=2.0;
    ParKeyVal *pkv2 = NULL;
    pkv2 = ParKeyVal_add(pkv2, "x", &xx, Free);
    pkv2 = ParKeyVal_add(pkv2, "y", &yy, Free);
    pkv2 = ParKeyVal_add(pkv2, "z", &zz, Free);
    unitTstResult("ParKeyVal", "OK");

    ParStore *ps = ParStore_new();
    assert(ParStore_nFixed(ps) == 0);
    assert(ParStore_nFree(ps) == 0);
    assert(ParStore_nGaussian(ps) == 0);

    double val, *ptr;
    ParamStatus pstat, pstat2;

    val = 12.3;
    ParStore_addFixedPar(ps, val, "x");
    ptr = ParStore_findPtr(ps, &pstat, "x");
    assert(*ptr == val);
	assert(pstat == Fixed);
    assert(ParStore_nFixed(ps) == 1);
    assert(ParStore_nFree(ps) == 0);
    assert(ParStore_nGaussian(ps) == 0);
    assert(ParStore_getFixed(ps, 0) == val);
    ParStore_sample(ps, ptr, val-1,val+1, rng);
    assert(*ptr == val);

    val = 23.4;
    ParStore_addFreePar(ps, val, 10.0, 30.0, "y");
    ptr = ParStore_findPtr(ps, &pstat, "y");
    assert(*ptr == val);
	assert(pstat == Free);
    assert(ParStore_nFixed(ps) == 1);
    assert(ParStore_nFree(ps) == 1);
    assert(ParStore_nGaussian(ps) == 0);
    assert(ParStore_getFree(ps, 0) == val);
    assert(ParStore_loFree(ps, 0) == 10.0);
    assert(ParStore_hiFree(ps, 0) == 30.0);
    ParStore_sample(ps, ptr, val-1,val+1, rng);
    assert(*ptr == val);

    val = 88.3;
    ParStore_addFixedPar(ps, val, "w");
    ptr = ParStore_findPtr(ps, &pstat, "w");
    assert(*ptr == val);
	assert(pstat == Fixed);
    assert(ParStore_nFixed(ps) == 2);
    assert(ParStore_nFree(ps) == 1);
    assert(ParStore_nGaussian(ps) == 0);
    assert(ParStore_getFixed(ps, 1) == val);
    ParStore_sample(ps, ptr, val-1,val+1, rng);
    assert(*ptr == val);

    val = -23.8;
    ParStore_addFreePar(ps, val, -100.0, 0.0, "z");
    ptr = ParStore_findPtr(ps, &pstat, "z");
    assert(*ptr == val);
	assert(pstat == Free);
    assert(ParStore_nFixed(ps) == 2);
    assert(ParStore_nFree(ps) == 2);
    assert(ParStore_nGaussian(ps) == 0);
    assert(ParStore_getFree(ps, 1) == val);
    assert(ParStore_loFree(ps, 1) == -100.0);
    assert(ParStore_hiFree(ps, 1) == 0.0);
    ParStore_sample(ps, ptr, val-1,val+1, rng);
    assert(*ptr == val);

    val = 0.8;
    ParStore_addFreePar(ps, val, 0.0, 1.0, "a");
    ptr = ParStore_findPtr(ps, &pstat, "a");
    assert(*ptr == val);
	assert(pstat == Free);
    assert(ParStore_nFixed(ps) == 2);
    assert(ParStore_nFree(ps) == 3);
    assert(ParStore_nGaussian(ps) == 0);
    assert(ParStore_getFree(ps, 2) == val);
    assert(ParStore_loFree(ps, 2) == 0.0);
    assert(ParStore_hiFree(ps, 2) == 1.0);
    ParStore_sample(ps, ptr, val-1,val+1, rng);
    assert(*ptr == val);

    double sd = 1.2;
    val = 2.3;
    ParStore_addGaussianPar(ps, val, sd, "gauss");
    ptr = ParStore_findPtr(ps, &pstat, "gauss");
    assert(*ptr = val);
	assert(pstat == Gaussian);
    assert(ParStore_nFixed(ps) == 2);
    assert(ParStore_nFree(ps) == 3);
    assert(ParStore_nGaussian(ps) == 1);
    assert(ParStore_getGaussian(ps, 0) == val);
    assert(*ptr == val);
    ParStore_sample(ps, ptr, val-1,val+1, rng);
    assert(ParStore_getGaussian(ps, 0) != val);
    assert(*ptr != val);

    strcpy(buff, "1 + 2*y + 1*a");
    assert(0 == ParStore_nConstrained(ps));
    ParStore_addConstrainedPar(ps, buff, "cnstr");
    assert(1 == ParStore_nConstrained(ps));
    ParStore_constrain(ps);

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
        size_t position1 = (size_t) ParStore_findPtr(ps, &pstat, name1);
        size_t position2 = (size_t) ParStore_findPtr(ps2,&pstat2, name1);
        assert(offset == position2 - position1);
    }

    for(i=0; i < ParStore_nFree(ps); ++i) {
        assert(ParStore_getFree(ps2, i) == ParStore_getFree(ps, i));
        assert(ParStore_loFree(ps2, i) == ParStore_loFree(ps, i));
        assert(ParStore_hiFree(ps2, i) == ParStore_hiFree(ps, i));
        const char *name1 = ParStore_getNameFree(ps, i);
        const char *name2 = ParStore_getNameFree(ps2, i);
        assert(0 == strcmp(name1, name2));
        size_t position1 = (size_t) ParStore_findPtr(ps, &pstat,name1);
        size_t position2 = (size_t) ParStore_findPtr(ps2,&pstat2,name1);
        assert(offset == position2 - position1);
    }
    if(verbose)
        ParStore_print(ps2, stdout);

    ParStore_free(ps);
    ParStore_free(ps2);
    ParKeyVal_free(pkv);
    ParKeyVal_free(pkv2);
    gsl_rng_free(rng);
    unitTstResult("ParStore", "OK");

    return 0;
}
