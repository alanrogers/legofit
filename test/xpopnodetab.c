/**
 * @file xpopnodetab.c
 * @author Alan R. Rogers
 * @brief Test popnodetab.c.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "popnodetab.h"
#include "misc.h"
#include "gptree.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {
    int verbose = 0;

	switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xpopnodetab [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xpopnodetab [-v]\n");
    }


    int nnodes = 3;
    PopNode pn[nnodes];
    const char *lbl[] = {"jack", "jill", "tom"};

    PopNodeTab *pnt = PopNodeTab_new();

    assert(NULL == PopNodeTab_get(pnt, "foo"));
    assert(0 == PopNodeTab_size(pnt));

    assert(0 == PopNodeTab_insert(pnt, lbl[0], &pn[0]));
    assert(1 == PopNodeTab_size(pnt));
    assert(&pn[0] == PopNodeTab_get(pnt, lbl[0]));
    assert(NULL == PopNodeTab_get(pnt, "foo"));

    assert(0 == PopNodeTab_insert(pnt, lbl[1], &pn[1]));
    assert(2 == PopNodeTab_size(pnt));
    assert(&pn[0] == PopNodeTab_get(pnt, lbl[0]));
    assert(&pn[1] == PopNodeTab_get(pnt, lbl[1]));
    assert(NULL == PopNodeTab_get(pnt, "foo"));

    assert(0 == PopNodeTab_insert(pnt, lbl[2], &pn[2]));
    assert(3 == PopNodeTab_size(pnt));
    assert(&pn[0] == PopNodeTab_get(pnt, lbl[0]));
    assert(&pn[1] == PopNodeTab_get(pnt, lbl[1]));
    assert(&pn[2] == PopNodeTab_get(pnt, lbl[2]));
    assert(NULL == PopNodeTab_get(pnt, "foo"));

    if(verbose)
        PopNodeTab_print(pnt);

    PopNodeTab_free(pnt);

    unitTstResult("PopNodeTab", "OK");

    return 0;
}
