/**
 * @file xidset.c
 * @author Alan R. Rogers
 * @brief Test idset.c.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "idset.h"
#include "misc.h"
#include "event.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif
 
int main(int argc, char **argv) {

    int         i, verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0) {
            fprintf(stderr, "usage: xidset [-v]\n");
            exit(1);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xidset [-v]\n");
        exit(1);
    }

    const int ntips = 4;

    IdSet *tip[ntips];

    for(i=0; i<ntips; ++i) {
        tip[i] = IdSet_newTip(1u << i);
        if(verbose) {
            IdSet_print(tip[i], stdout);
            putchar('\n');
        }
        IdSet_sanityCheck(tip[i], __FILE__, __LINE__);
    }

    int nIds = 5;
    tipId_t tid[5] = {1, 2, 4, 8, 16};
    double pr = 0.5;

    IdSet *a = IdSet_new(nIds, tid, pr);
    IdSet_sanityCheck(a, __FILE__, __LINE__);

    if(verbose) {
        printf("A set with 5 ids\n");
        IdSet_print(a, stdout);
        putchar('\n');
    }

    IdSet *b = IdSet_dup(a);
    IdSet_sanityCheck(b, __FILE__, __LINE__);

    assert(0 == IdSet_cmp(a, b));

    if(verbose) {
        printf("Two duplicate IdSet objects\n");
        printf("***a:\n");
        IdSet_print(a, stdout);
        putchar('\n');
        printf("***b:\n");
        IdSet_print(b, stdout);
        putchar('\n');
    }

    nIds = 3;
    tipId_t tid2[3] = {32,64,128};
    pr = 0.5;
    IdSet *c = IdSet_new(nIds, tid2, pr);
    IdSet_sanityCheck(c, __FILE__, __LINE__);

    assert(0 < IdSet_cmp(a, c));

    IdSet_addMigEvent(a, 1,1, 0.1);
    IdSet_addMigEvent(a, 2,2, 0.2);
    IdSet_addMigEvent(c, 3,3, 0.3);

    IdSet *d = IdSet_join(a, c, 0, NULL);
    assert(d);
    IdSet_sanityCheck(d, __FILE__, __LINE__);
    if(verbose) {
        printf("a join c = d\n");
        printf("***a:\n");
        IdSet_print(a, stdout);
        putchar('\n');
        printf("***c:\n");
        IdSet_print(c, stdout);
        putchar('\n');
        printf("***d:\n");
        IdSet_print(d, stdout);
        putchar('\n');
    }

    // construct "e", which is mutually exclusive with "a".
    IdSet *e = IdSet_new(nIds, tid2, pr);
    IdSet_sanityCheck(e, __FILE__, __LINE__);
    IdSet_addMigEvent(e, 1, 2, 0.1);

    // can't join IdSet objects with mutually exclusive
    // migration histories.
    assert(NULL == IdSet_join(a, e, 0, NULL));
    
    for(i=0; i<ntips; ++i)
        IdSet_free(tip[i]);
    IdSet_free(a);
    IdSet_free(b);
    IdSet_free(c);
    IdSet_free(d);
    IdSet_free(e);

    unitTstResult("IdSet", "OK");
    
    return 0;
}
