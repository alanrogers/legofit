/**
 * @file xu64u64map.c
 * @author Alan R. Rogers
 * @brief Test u64u64.c.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "u64u64map.h"
#include "misc.h"
#include "error.h"
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
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
        if(strncmp(argv[1], "-v", 2) != 0) {
            fprintf(stderr, "usage: xu64u64map [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xu64u64map [-v]\n");
        exit(EXIT_FAILURE);
    }


    U64U64Map *map = U64U64Map_new(8);

    const int nvals = 50;
    uint64_t key[nvals], value[nvals], val;
    int status;

    for(int i=0; i < nvals; ++i) {
        key[i] = rand();
        value[i] = rand();
        status = U64U64Map_insert(map, key[i], value[i]);
        assert(status == 0);
    }

    unsigned long size = U64U64Map_size(map);
    assert(nvals == size);

    uint64_t kk[size];
    status = U64U64Map_keys(map, size, kk);
    assert(status == 0);

    if(verbose) {
        for(int i=0; i<size; ++i) {
            val = U64U64Map_get(map, kk[i], &status);
            assert(status==0);
            printf("%llu => %llu\n", kk[i], val);
        }                 
    }

    status = U64U64Map_keys(map, size/2, kk);
    assert(status == BUFFER_OVERFLOW);

    for(int i=0; i < nvals; ++i) {
        val = U64U64Map_get(map, key[i], &status);
        assert(status==0);
        assert(val == value[i]);
    }

    // find a k that is not a key in the map
    uint64_t k, k_in_map;
    do{
        k = rand();
        k_in_map = 0;
        for(int i=0; i<nvals; ++i) {
            if(k == key[i]) {
                k_in_map = 1;
                break;
            }
        }
    }while(k_in_map);

    val = U64U64Map_get(map, k, &status);
    assert(status == 1);

    U64U64Map_free(map);

    unitTstResult("U64U64Map", "OK");

    return 0;
}
