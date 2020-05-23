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
#include "network.h"
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

    int verbose = 0;

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

    Network_init(SIM);

    Param par;
    Param_init(&par, "name", 123.4, 100.0, 200.0, FREE|TWON);
    assert(strcmp(par.name, "name") == 0);
    assert(par.value == 123.4);
    assert(par.low == 100.0);
    assert(par.high == 200.0);

    if(verbose)
        Param_print(&par, stdout);

    Param_sanityCheck(&par, __FILE__, __LINE__);

    if(verbose)
        Param_print(&par, stdout);

    Param_freePtrs(&par);

    unitTstResult("Param", "OK");

    return 0;

 usage:
    fprintf(stderr, "usage xparam [-v]\n");
    exit(EXIT_FAILURE);
}
