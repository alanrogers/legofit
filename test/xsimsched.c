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

int main(void) {

    int itr[] = {100, 20, 300};
    int reps[] = {1000, 2000, 3000};

    SimSched *ss = SimSched_new();
    assert(SimSched_empty(ss));
    SimSched_append(ss, itr[0], reps[0]);
    assert(!SimSched_empty(ss));
    SimSched_append(ss, itr[1], reps[1]);
    assert(!SimSched_empty(ss));
    SimSched_append(ss, itr[2], reps[2]);
    assert(!SimSched_empty(ss));

    assert(SimSched_getOptItr(ss) == 100);
    assert(SimSched_getSimReps(ss) == 1000);
    SimSched_next(ss);
    assert(SimSched_getOptItr(ss) == 20);
    assert(SimSched_getSimReps(ss) == 2000);
    SimSched_next(ss);
    assert(SimSched_getOptItr(ss) == 300);
    assert(SimSched_getSimReps(ss) == 3000);
    SimSched_next(ss);
    assert(SimSched_empty(ss));
    SimSched_free(ss);

    unitTstResult("SimSched", "OK");

    return 0;
}
