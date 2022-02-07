/**
 * @file xidsetset.c
 * @author Alan R. Rogers
 * @brief Test idset_tbl.c.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "idsetset.h"
#include "idset.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif
 
int main(int argc, char **argv) {

    int i, verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0) {
            fprintf(stderr, "usage: xidset_tbl [-v]\n");
            exit(1);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xidset_tbl [-v]\n");
        exit(1);
    }

    IdSetSet *iss = IdSetSet_new(0); // will round up to next pwr of 2
    CHECKMEM(iss);
    IdSetSet_sanityCheck(iss, __FILE__, __LINE__);

    int status;
    const int ntips = 4;

    IdSet *tip[ntips];

    for(i=0; i<ntips; ++i) {
        tip[i] = IdSet_newTip(1u << i);
        IdSet_sanityCheck(tip[i], __FILE__, __LINE__);
        status = IdSetSet_add(iss, tip[i]);
        if(status)
            ERR(status, "bad return from IdSetSet_add");
        IdSetSet_sanityCheck(iss, __FILE__, __LINE__);
    }

    int nIds = 5;
    tipId_t tid[5] = {1, 2, 4, 8, 16};
    double pr = 0.5;
    long unsigned event = 0;
    long unsigned outcome = 0;

    IdSet *a = IdSet_new(nIds, tid, EventLst_insert(NULL,event,outcome,pr));
    IdSet_sanityCheck(a, __FILE__, __LINE__);

    status = IdSetSet_add(iss, a);
    if(status)
        ERR(status, "bad return from IdSetSet_add");
    IdSetSet_sanityCheck(iss, __FILE__, __LINE__);

    IdSet *b = IdSet_dup(a);
    IdSet_sanityCheck(b, __FILE__, __LINE__);

    assert(0 == IdSet_cmp(a, b));

    nIds = 3;
    tipId_t tid2[3] = {32,64,128};
    pr = 0.5;
    IdSet *c = IdSet_new(nIds, tid2, EventLst_insert(NULL,event,outcome,pr));
    IdSet_sanityCheck(c, __FILE__, __LINE__);

    assert(0 < IdSet_cmp(a, c));

    IdSet_addEvent(a, 1,1, 0.1);
    IdSet_addEvent(a, 2,2, 0.2);
    IdSet_addEvent(c, 3,3, 0.3);

    status = IdSetSet_add(iss, c);
    if(status)
        ERR(status, "bad return from IdSetSet_add");

    IdSet *d = IdSet_join(a, c, 0, NULL);
    assert(d);
    IdSet_sanityCheck(d, __FILE__, __LINE__);

    status = IdSetSet_add(iss, d);
    if(status)
        ERR(status, "bad return from IdSetSet_add");

    // construct "e", which is mutually exclusive with "a".
    IdSet *e = IdSet_new(nIds, tid2, EventLst_insert(NULL,event,outcome,pr));
    IdSet_sanityCheck(e, __FILE__, __LINE__);
    IdSet_addEvent(e, 1, 2, 0.1);

    status = IdSetSet_add(iss, e);
    if(status)
        ERR(status, "bad return from IdSetSet_add");

    // can't join IdSet objects with mutually exclusive
    // migration histories.
    assert(NULL == IdSet_join(a, e, 0, NULL));

    IdSetSet_rewind(iss);
    IdSetSet_sanityCheck(iss, __FILE__, __LINE__);
    for(a=IdSetSet_next(iss); a; a=IdSetSet_next(iss))
        IdSet_sanityCheck(a, __FILE__,__LINE__);

    if(verbose) {
        IdSetSet_rewind(iss);
        for(a=IdSetSet_next(iss); a; a=IdSetSet_next(iss))
            IdSet_print(a, stderr);
        putc('\n', stderr);
    }

    IdSet_free(b);

    IdSetSet_sanityCheck(iss, __FILE__, __LINE__);
    IdSetSet_free_deep(iss);

    unitTstResult("IdSetSet", "OK");
    
    return 0;
}
