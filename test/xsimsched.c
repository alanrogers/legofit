/**
 * @file xsimsched.c
 * @author Alan R. Rogers
 * @brief Test parstore.c.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "simsched.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

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
            eprintf("usage: xsimsched [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xsimsched [-v]\n");
    }

    long itr[] = {100L, 20L, 300L};
    long reps[] = {1000L, 2000L, 3000L};

    SimSched *ss = SimSched_new();
    assert(SimSched_empty(ss));
    SimSched_append(ss, itr[0], reps[0]);
    assert(!SimSched_empty(ss));
    SimSched_append(ss, itr[1], reps[1]);
    assert(!SimSched_empty(ss));
    SimSched_append(ss, itr[2], reps[2]);
    assert(!SimSched_empty(ss));

    if(verbose)
        SimSched_print(ss, stdout);

    SimSched *ss2 = SimSched_dup(ss);

    assert(SimSched_getOptItr(ss) == 100L);
    assert(SimSched_getOptItr(ss2) == 100L);
    assert(SimSched_getSimReps(ss) == 1000L);
    assert(SimSched_getSimReps(ss2) == 1000L);
    SimSched_next(ss);
    SimSched_next(ss2);
    assert(SimSched_getOptItr(ss) == 20L);
    assert(SimSched_getOptItr(ss2) == 20L);
    assert(SimSched_getSimReps(ss) == 2000L);
    assert(SimSched_getSimReps(ss2) == 2000L);
    SimSched_next(ss);
    SimSched_next(ss2);
    assert(SimSched_getOptItr(ss) == 300L);
    assert(SimSched_getOptItr(ss2) == 300L);
    assert(SimSched_getSimReps(ss) == 3000L);
    assert(SimSched_getSimReps(ss2) == 3000L);
    SimSched_next(ss);
    SimSched_next(ss2);
    assert(SimSched_empty(ss));
    assert(SimSched_empty(ss2));
    SimSched_free(ss);
    SimSched_free(ss2);


    unitTstResult("SimSched", "OK");

    return 0;
}
