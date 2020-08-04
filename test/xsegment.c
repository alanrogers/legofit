/**
 * @file xsegment.c
 * @author Alan R. Rogers
 * @brief Test segment.c.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "branchtab.h"
#include "error.h"
#include "misc.h"
#include "param.h"
#include "parstore.h"
#include "ptrptrmap.h"
#include "ptrqueue.h"
#include "segment.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char *argv[]) {

    int status, verbose = 0;

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


    PtrQueue *fixedQ = PtrQueue_new();
    PtrQueue *freeQ = PtrQueue_new();
    PtrQueue *constrQ = PtrQueue_new();

    Param *par;

    par = Param_new("zero", 0.0, 0.0, 0.0, TIME|FIXED, NULL);
    PtrQueue_push(fixedQ, par);

    par = Param_new("one", 1.0, 1.0, 1.0, TWON|FIXED, NULL);
    PtrQueue_push(fixedQ, par);

    par = Param_new("Nab", 3.0, 0.0, 100.0, TWON|FREE, NULL);
    PtrQueue_push(freeQ, par);

    par = Param_new("Tab", 2.0, 0.0, 100.0, TIME|FREE, NULL);
    PtrQueue_push(freeQ, par);

    par = Param_new("Tmig", 1.0, 1.0, 1.0, TIME|FIXED, NULL);
    PtrQueue_push(fixedQ, par);

    par = Param_new("mix", 0.02, 0.02, 0.02, MIXFRAC|FIXED, NULL);
    PtrQueue_push(fixedQ, par);

    par = Param_new("Nabc", 3.0, 0.0, 100.0, TWON|FREE, NULL);
    PtrQueue_push(freeQ, par);

    par = Param_new("Tabc", 4.0, -DBL_MAX, DBL_MAX, TIME|CONSTRAINED,
                    "Tab + Nab*Nabc");
    PtrQueue_push(constrQ, par);

    ParStore *ps = ParStore_new(fixedQ, freeQ, constrQ);

    if(verbose)
        ParStore_print(ps, stderr);

    Bounds bnd = {
        .lo_twoN = 0.0,
        .hi_twoN = 1e12,
        .lo_t = 0,
        .hi_t = 1e10
    };

    Segment *a, *b, *b2, *c, *c2, *ab, *abc;
    int ni, ti, mi;
#if 0    
    tipId_t     ida = 0;
    tipId_t     idb = 1;
    tipId_t     idc = 2;
#endif    

    ni = ParStore_getIndex(ps, "one");
    assert(ni >= 0);
    ti = ParStore_getIndex(ps, "zero");
    assert(ti >= 0);

    a = Segment_new(ni, ti, ps);
    assert(a);

    b = Segment_new(ni, ti, ps);
    assert(b);

    c = Segment_new(ni, ti, ps);
    assert(c);
    
    Segment_newSample(a, 0);
    Segment_newSample(b, 1);
    Segment_newSample(c, 2);

    ti = ParStore_getIndex(ps, "Tmig");
    assert(ti >= 0);
    b2 = Segment_new(ni, ti, ps);
    c2 = Segment_new(ni, ti, ps);
    assert(b2);
    assert(c2);

    ni = ParStore_getIndex(ps, "Nab");
    assert(ni >= 0);
    ti = ParStore_getIndex(ps, "Tab");
    assert(ti >= 0);
    ab = Segment_new(ni, ti, ps);
    assert(ab);

    ni = ParStore_getIndex(ps, "Nabc");
    assert(ni >= 0);
    ti = ParStore_getIndex(ps, "Tabc");
    assert(ti >= 0);
    abc = Segment_new(ni, ti, ps);
    assert(abc);
    
    status = Segment_addChild(ab, a);
    assert(status == 0);

    mi = ParStore_getIndex(ps, "mix");
    assert(mi >= 0);
    status = Segment_mix(b, mi, c2, b2, ps);

    status = Segment_addChild(c2, c);
    assert(status == 0);

    status = Segment_addChild(ab, b2);
    assert(status == 0);
    
    status = Segment_addChild(abc, ab);
    assert(status == 0);

    status = Segment_addChild(abc, c2);
    assert(status == 0);

    assert(abc == Segment_root(a));
    assert(abc == Segment_root(b));
    assert(abc == Segment_root(b2));
    assert(abc == Segment_root(c));
    assert(abc == Segment_root(c2));
    assert(abc == Segment_root(ab));
    assert(abc == Segment_root(abc));

    assert(Segment_feasible(abc, bnd, verbose));

    PtrPtrMap *ppm = PtrPtrMap_new();
    Segment *duproot = Segment_dup(abc, ppm);
    CHECKMEM(duproot);
    assert(Segment_feasible(duproot, bnd, verbose));
    PtrPtrMap_free(ppm);

    assert(Segment_equals(abc, duproot));

    if(verbose)
        Segment_print(stdout, abc, 0);

    BranchTab *bt = BranchTab_new();

    Segment_unvisit(abc);
    status = Segment_coalesce(abc, 1, bt);
    if(status) {
        fprintf(stderr, "%s:%d: Segment_coalesce returned %d\n",
                __FILE__,__LINE__,status);
        exit(EXIT_FAILURE);
    }
    BranchTab_print(bt, stdout);

    Segment_free(abc);
    Segment_free(duproot);
    
    unitTstResult("Segment", "OK");

    BranchTab_free(bt);

    return 0;

 usage:
    fprintf(stderr, "usage xsegment [-v]\n");
    exit(EXIT_FAILURE);
}
