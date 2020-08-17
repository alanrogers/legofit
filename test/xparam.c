/**
 * @file xparam.c
 * @author Alan R. Rogers
 * @brief Test param.c.
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "param.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <limits.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char *argv[]) {

    int verbose = 0, status;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            goto usage;
        verbose = 1;
        break;
    default:
        goto usage;
    }

    time_t      currtime = time(NULL);
    unsigned    seed = currtime % UINT_MAX;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    CHECKMEM(rng);
    gsl_rng_set(rng, seed);

    // par lives in [100, 200]
    Param *par = Param_new("name", 123.4, 100.0, 200.0, FREE|TWON, NULL);
    assert(strcmp(par->name, "name") == 0);
    assert(par->value == 123.4);
    assert(par->low == 100.0);
    assert(par->high == 200.0);
    assert(Param_isFree(par));

    double v;
    for(int i=0; i < 4; ++i) {
        v = Param_getTrialValue(par, rng);
        assert(v >= 100.0);
        assert(v < 200.0);
        assert(v != 100.0);
    }

    // par2 lives in [0, 1]
    Param *par2 = Param_new("foo", 0.5, 0.0, 1.0, FIXED|MIXFRAC, NULL);
    Param_sanityCheck(par2, __FILE__, __LINE__);
    assert(!Param_isFree(par2));

    assert(0 != Param_compare(par, par2));

    // try assigning a value out of range
    assert(EDOM == Param_setValue(par2, 99.99));

    // try assigning a value in range
    status = Param_setValue(par2, 0.25);
    assert(status != EDOM);
    assert(0 == status);
    assert(0.25 == Param_getValue(par2));

    // make par2 equal to par
    Param_copy(par2, par);
    assert(0 == Param_compare(par, par2));
    assert(Param_isFree(par2));

    Param_freePtrs(par);
    free(par);

    StrPtrMap *pars = StrPtrMap_new();

    char formula[100];
    double w=1.0, x=2.0, y=3.0, z=4.0;
    StrPtrMap_insert(pars, "w", &w);
    StrPtrMap_insert(pars, "x", &x);
    StrPtrMap_insert(pars, "y", &y);
    StrPtrMap_insert(pars, "z", &z);
    sprintf(formula, "%s", "w + x*y - z");

    // A constrained Param
    par = Param_new("constrained", 4.0, -DBL_MAX, DBL_MAX,
                    CONSTRAINED|TWON, formula);
    assert(!Param_isFree(par));
    assert(4.0 == Param_getValue(par));
    Param_compileConstraint(par, pars);
    Param_constrain(par);
    assert(3.0 == Param_getValue(par));

    // Change one of the underlying variables
    w += 1;
    Param_constrain(par);
    assert(4.0 == Param_getValue(par));

    Param_move(par2, par);
    free(par);
    assert(par2->type == (CONSTRAINED|TWON));
    assert(4.0 == Param_getValue(par2));

    // Change one of the underlying variables
    w += 1;
    Param_constrain(par2);
    assert(5.0 == Param_getValue(par2));
    Param_freePtrs(par2);
    free(par2);

    unitTstResult("Param", "OK");

    gsl_rng_free(rng);
    te_free_variables(pars);

    return 0;

 usage:
    fprintf(stderr, "usage xparam [-v]\n");
    exit(EXIT_FAILURE);
}
