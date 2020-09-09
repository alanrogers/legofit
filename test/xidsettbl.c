/**
 * @file xidsettbl.c
 * @author Alan R. Rogers
 * @brief Test idset_tbl.c.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "idsettbl.h"
#include "idset.h"
#include "misc.h"
#include "migoutcome.h"
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

    IdSetTbl *ist = IdSetTbl_new(0); // will round up to next pwr of 2
    CHECKMEM(ist);

    int status;
    const int ntips = 4;

    IdSet *tip[ntips];

    for(i=0; i<ntips; ++i) {
        tip[i] = IdSet_newTip(1u << i);
        IdSet_sanityCheck(tip[i], __FILE__, __LINE__);
        status = IdSetTbl_add(ist, tip[i]);
        switch(status) {
        case ENOMEM:
            fprintf(stderr,"%s:%d: bad alloc\n",__FILE__,__LINE__);
            exit(EXIT_FAILURE);
        case 0:
            break;
        default:
            fprintf(stderr,"%s:%d: unknown error\n",__FILE__,__LINE__);
            exit(EXIT_FAILURE);
        }
    }

    int nIds = 5;
    tipId_t tid[5] = {1, 2, 4, 8, 16};
    double pr = 0.5;

    IdSet *a = IdSet_new(nIds, tid, pr);
    IdSet_sanityCheck(a, __FILE__, __LINE__);

    status = IdSetTbl_add(ist, a);
    assert(status == 0);

    IdSet *b = IdSet_dup(a);
    IdSet_sanityCheck(b, __FILE__, __LINE__);

    assert(0 == IdSet_cmp(a, b));

    status = IdSetTbl_add(ist, b);
    assert(status == 0);

    nIds = 3;
    tipId_t tid2[3] = {32,64,128};
    pr = 0.5;
    IdSet *c = IdSet_new(nIds, tid2, pr);
    IdSet_sanityCheck(c, __FILE__, __LINE__);

    assert(0 < IdSet_cmp(a, c));

    IdSet_addMigEvent(a, 1,1, 0.1);
    IdSet_addMigEvent(a, 2,2, 0.2);
    IdSet_addMigEvent(c, 3,3, 0.3);

    status = IdSetTbl_add(ist, c);
    assert(status == 0);

    IdSet *d = IdSet_join(a, c, 0, NULL);
    assert(d);
    IdSet_sanityCheck(d, __FILE__, __LINE__);

    status = IdSetTbl_add(ist, d);
    assert(status == 0);

    // construct "e", which is mutually exclusive with "a".
    IdSet *e = IdSet_new(nIds, tid2, pr);
    IdSet_sanityCheck(e, __FILE__, __LINE__);
    IdSet_addMigEvent(e, 1, 2, 0.1);

    status = IdSetTbl_add(ist, e);
    assert(status == 0);

    // can't join IdSet objects with mutually exclusive
    // migration histories.
    assert(NULL == IdSet_join(a, e, 0, NULL));

    IdSetTbl_rewind(ist);
    for(a=IdSetTbl_next(ist); a; a=IdSetTbl_next(ist))
        IdSet_sanityCheck(a, __FILE__,__LINE__);

    if(verbose) {
        IdSetTbl_rewind(ist);
        for(a=IdSetTbl_next(ist); a; a=IdSetTbl_next(ist))
            IdSet_print(a, stderr);
        putc('\n', stderr);
    }

    IdSetTbl_rewind(ist);
    for(a=IdSetTbl_next(ist); a; a=IdSetTbl_next(ist))
        IdSet_free(a);
    

    IdSetTbl_free(ist);

    unitTstResult("IdSetTbl", "OK");
    
    return 0;
}
