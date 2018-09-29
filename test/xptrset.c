/**
 * @file xptrset.c
 * @author Alan R. Rogers
 * @brief Unit test for PtrSet
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "ptrset.h"
#include "misc.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>

#ifdef NDEBUG
#  error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {
    int         verbose = 0;

    if(argc==2 && 0==strcmp("-v", argv[1]))
        verbose = 1;
    else if(argc != 1) {
        fprintf(stderr,"usage: xdtnorm [-v]\n");
        exit(1);
    }

    double x=1.0, y=2.0, z=3.0;

    PtrSet *ps = PtrSet_new();

    assert( !PtrSet_exists(ps, &x) );

    PtrSet_insert(ps, &y);
    assert( !PtrSet_exists(ps, &x) );
    assert( PtrSet_exists(ps, &y) );
    assert( !PtrSet_exists(ps, &z) );

    PtrSet_insert(ps, &x);
    assert( PtrSet_exists(ps, &x) );
    assert( PtrSet_exists(ps, &y) );
    assert( !PtrSet_exists(ps, &z) );

    PtrSet_insert(ps, &z);
    assert( PtrSet_exists(ps, &x) );
    assert( PtrSet_exists(ps, &y) );
    assert( PtrSet_exists(ps, &z) );

    // Inserting the same pointer multiple times should have no effect
    PtrSet_insert(ps, &x);
    PtrSet_insert(ps, &y);
    PtrSet_insert(ps, &z);

    if(verbose) {
        printf("PtrSet list should have 3 items:");
        PtrSet_print(ps, stdout);
    }

    unitTstResult("PtrSet", "OK");
    return 0;
}
