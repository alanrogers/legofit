/**
 * @file xptrlst.c
 * @author Alan R. Rogers
 * @brief Test ptrlst.c.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

#include "ptrlst.h"
#include "misc.h"
#include <stdio.h>
#include <assert.h>
 
int main(void) {

    int a=1, b=2, c=3, d=4, status;
    PtrLst *v = PtrLst_new();
    assert(0 == PtrLst_length(v));

    status = PtrLst_push(v, &a);
    assert(status == 0);
    assert(1 == PtrLst_length(v));

    status = PtrLst_push(v, &b);
    assert(status == 0);
    assert(2 == PtrLst_length(v));
    
    status = PtrLst_push(v, &c);
    assert(status == 0);
    assert(3 == PtrLst_length(v));
    
    status = PtrLst_push(v, &d);
    assert(status == 0);
    assert(4 == PtrLst_length(v));

    int *p = PtrLst_pop(v);
    assert(p == &d);
    assert(3 == PtrLst_length(v));

    p = PtrLst_pop(v);
    assert(p == &c);
    assert(2 == PtrLst_length(v));

    p = PtrLst_pop(v);
    assert(p == &b);
    assert(1 == PtrLst_length(v));

    p = PtrLst_pop(v);
    assert(p == &a);
    assert(0 == PtrLst_length(v));
    
    p = PtrLst_pop(v);
    assert(p == NULL);

    PtrLst_free(v);

    unitTstResult("PtrLst", "OK");
    return 0;
}
