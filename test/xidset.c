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
        if(i>0)
            assert(0 != IdSet_cmp(tip[i-1], tip[i]));
        if(verbose)
            IdSet_print(tip[i], stdout);
        IdSet_sanityCheck(tip[i], __FILE__, __LINE__);
    }

    int nIds = 5;
    tipId_t tid[5] = {1, 2, 3, 4, 5};
    double pr = 0.5;

    IdSet *a = IdSet_new(nIds, tid, pr);
    IdSet_sanityCheck(a, __FILE__, __LINE__);

    if(verbose) {
        printf("A set with 5 ids\n");
        IdSet_print(a, stdout);
    }

    IdSet *b = IdSet_dup(a);
    IdSet_sanityCheck(b, __FILE__, __LINE__);

    assert(0 == IdSet_cmp(a, b));

    for(i=0; i<ntips; ++i)
        IdSet_free(tip[i]);
    IdSet_free(a);
    IdSet_free(b);

    unitTstResult("IdSet", "OK");
    
    return 0;
}
