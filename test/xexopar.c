/**
 * @file xexopar.c
 * @author Alan R. Rogers
 * @brief Test tokenizer.c.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "exopar.h"
#include "misc.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <limits.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {

    int         verbose = 0;
    time_t      currtime = time(NULL);
    unsigned    baseSeed = currtime % UINT_MAX;
    gsl_rng    *rng = gsl_rng_alloc(gsl_rng_taus);
    CHECKMEM(rng);
    gsl_rng_set(rng, baseSeed);

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xexopar [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xexopar [-v]\n");
    }

    double x=0, y=0, z=0;
    double x2=123.45, y2=456.78, z2=910.11;

    if(verbose)
        printf("init: x=%lf y=%lf z=%lf\n", x, y, z);

    ExoPar *ep = ExoPar_new();
    CHECKMEM(ep);

    ExoPar_add(ep, &x, 1.0, 1.0);
    ExoPar_add(ep, &y, 2.0, 1.0);
    ExoPar_add(ep, &z, 3.0, 1.0);
    ExoPar_freeze(ep);

    assert(x == 1.0);
    assert(y == 2.0);
    assert(z == 3.0);
    if(verbose)
        printf("after ExoPar_add: x=%lf y=%lf z=%lf\n",
               x, y, z);

    assert(1 == ExoPar_sample(ep, &x2, 1.1, 1.2, rng));

    assert(0 == ExoPar_sample(ep, &x, -10.0, 10.0, rng));
    assert(0 == ExoPar_sample(ep, &y, -10.0, 10.0, rng));
    assert(0 == ExoPar_sample(ep, &z, -10.0, 10.0, rng));
    assert(-10.0 <= x && x <= 10.0);
    assert(-10.0 <= y && y <= 10.0);
    assert(-10.0 <= z && z <= 10.0);
    if(verbose)
        printf("after ExoPar_sample: x=%lf y=%lf z=%lf\n",
               x, y, z);

    assert(0 == ExoPar_sample(ep, &x, 1.1, 1.2, rng));
    assert(0 == ExoPar_sample(ep, &y, 1.1, 1.2, rng));
    assert(0 == ExoPar_sample(ep, &z, 1.1, 1.2, rng));
    assert(1.1 <= x && x <= 1.2);
    assert(1.1 <= y && y <= 1.2);
    assert(1.1 <= z && z <= 1.2);
    if(verbose)
        printf("after ExoPar_sample: x=%lf y=%lf z=%lf\n",
               x, y, z);

    size_t offset;
    int sign;
    if(&x2 > &x) {
        offset = ((size_t) &x2) - ((size_t) &x);
        sign = 1;
    }else{
        offset = ((size_t) &x) - ((size_t) &x2);
        sign = -1;
    }
    ExoPar_shiftPtrs(ep, offset, sign);

    assert(1 == ExoPar_sample(ep, &x, 1.1, 1.2, rng));

    assert(0 == ExoPar_sample(ep, &x2, 1.1, 1.2, rng));
    assert(0 == ExoPar_sample(ep, &y2, 1.1, 1.2, rng));
    assert(0 == ExoPar_sample(ep, &z2, 1.1, 1.2, rng));
    assert(1.1 <= x2 && x2 <= 1.2);
    assert(1.1 <= y2 && y2 <= 1.2);
    assert(1.1 <= z2 && z2 <= 1.2);
    if(verbose)
        printf("after shiftPtrs: x2=%lf y2=%lf z2=%lf\n",
               x2, y2, z2);

    ExoPar_free(ep);
    gsl_rng_free(rng);
    unitTstResult("ExoPar", "OK");
    return 0;
}
