/**
 * @file xptru32map.c
 * @author Alan R. Rogers
 * @brief Test ptru32map.c.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "ptru32map.h"
#include "misc.h"
#include "error.h"
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {
    int verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0) {
            fprintf(stderr, "usage: xptru32map [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xptru32map [-v]\n");
        exit(EXIT_FAILURE);
    }


    PtrU32Map *map = PtrU32Map_new();

    const int nvals = 50;
    unsigned key[nvals] = {0};
    uint32_t  value[nvals], val;
    int status;

    for(int i=0; i < nvals; ++i) {
        value[i] = rand();
        status = PtrU32Map_insert(map, key+i, value[i]);
        assert(status == 0);
    }

    unsigned long size = PtrU32Map_size(map);
    assert(nvals == size);

    void *kk[size];
    status = PtrU32Map_keys(map, size, kk);
    assert(status == 0);

    if(verbose) {
        for(int i=0; i<size; ++i) {
            val = PtrU32Map_get(map, kk[i], &status);
            assert(status==0);
            printf("%p => %u\n", kk[i], val);
        }                 
    }

    status = PtrU32Map_keys(map, size/2, kk);
    assert(status == BUFFER_OVERFLOW);

    for(int i=0; i < nvals; ++i) {
        val = PtrU32Map_get(map, key+i, &status);
        assert(status==0);
        assert(val == value[i]);
    }

    // k is not a key in the map
    unsigned *k = key + nvals;

    val = PtrU32Map_get(map, k, &status);
    assert(status == 1);

    PtrU32Map_free(map);

    unitTstResult("PtrU32Map", "OK");

    return 0;
}
