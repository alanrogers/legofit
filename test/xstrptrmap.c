/**
 * @file xstrptrmap.c
 * @author Alan R. Rogers
 * @brief Test strptrmap.c.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "strptrmap.h"
#include "misc.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int vptrcmp(const void *va, const void *vb);

// compare void pointers
int vptrcmp(const void *va, const void *vb) {
    const int * const * a = (const int * const *) va;
    const int * const * b = (const int * const *) vb;
    if(*a > *b)
        return 1;
    if(*a < *b)
        return -1;
    return 0;
}

int main(int argc, char **argv) {
    int verbose = 0;

	switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xstrptrmap [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xstrptrmap [-v]\n");
    }


    int val[3];
    const char *lbl[] = {"jack", "jill", "tom"};

    StrPtrMap *spm = StrPtrMap_new();

    assert(NULL == StrPtrMap_get(spm, "foo"));
    assert(0 == StrPtrMap_size(spm));

    assert(0 == StrPtrMap_insert(spm, lbl[0], &val[0]));
    assert(1 == StrPtrMap_size(spm));
    assert(&val[0] == StrPtrMap_get(spm, lbl[0]));
    assert(NULL == StrPtrMap_get(spm, "foo"));

    assert(0 == StrPtrMap_insert(spm, lbl[1], &val[1]));
    assert(2 == StrPtrMap_size(spm));
    assert(&val[0] == StrPtrMap_get(spm, lbl[0]));
    assert(&val[1] == StrPtrMap_get(spm, lbl[1]));
    assert(NULL == StrPtrMap_get(spm, "foo"));

    assert(0 == StrPtrMap_insert(spm, lbl[2], &val[2]));
    assert(3 == StrPtrMap_size(spm));
    assert(&val[0] == StrPtrMap_get(spm, lbl[0]));
    assert(&val[1] == StrPtrMap_get(spm, lbl[1]));
    assert(&val[2] == StrPtrMap_get(spm, lbl[2]));
    assert(NULL == StrPtrMap_get(spm, "foo"));

    void *ptr[3];
    StrPtrMap_ptrArray(spm, 3, ptr);

    // ptr contains the pointers in hashmap order. Use qsort
    // to put them into increasing order.
    qsort(ptr, 3, sizeof(void *), vptrcmp);

    assert(ptr[0] == (void *)(&val[0]));
    assert(ptr[1] == (void *)(&val[1]));
    assert(ptr[2] == (void *)(&val[2]));

    if(verbose)
        StrPtrMap_print(spm);

    StrPtrMap_free(spm);

    unitTstResult("StrPtrMap", "OK");

    return 0;
}
