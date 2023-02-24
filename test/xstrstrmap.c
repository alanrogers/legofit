/**
 * @file xstrstrmap.c
 * @author Alan R. Rogers
 * @brief Test strstrmap.c.
 * @copyright Copyright (c) 2023, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "strstrmap.h"
#include "misc.h"
#include "error.h"
#include <stdio.h>
#include <assert.h>
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
            fprintf(stderr, "usage: xstrstrmap [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xstrstrmap [-v]\n");
        exit(EXIT_FAILURE);
    }


    // Allocate small hash map so that it will be resized.
    StrStrMap *map = StrStrMap_new(8);

    assert(0 == StrStrMap_hasKey(map, "asdf"));

    const int nvals = 25;
    char *key[nvals];
    const char *val;
    int status;

    for(int i=0; i < nvals; ++i) {
        key[i] = malloc(2);
        key[i][0] = 'a' + i;
        key[i][1] = '\0';
        status = StrStrMap_insert(map, key[i], key[i]);
        assert(status == 0);
    }

    for(int i=0; i < nvals; ++i) {
        assert(1 == StrStrMap_hasKey(map, key[i]));
    }

    const char *longkey = "asdfasfadfasf asdfasfda fasdf adf asdf asdf ";
    status = StrStrMap_insert(map, longkey, "longkey");
    assert(status == 0);
    assert(1 == StrStrMap_hasKey(map, longkey));
    
    unsigned long size = StrStrMap_size(map);
    assert(nvals+1 == size);

    char * kk[size];
    status = StrStrMap_keys(map, size, kk);
    assert(status == 0);

    if(verbose) {
        for(int i=0; i<size; ++i) {
            val = StrStrMap_get(map, kk[i], &status);
            assert(status==0);
            printf("%s => %s\n", kk[i], val);
        }                 
    }
    for(int i=0; i < size; ++i)
        free(kk[i]);

    char * kkk[size/2];
    status = StrStrMap_keys(map, size/2, kkk);
    assert(status == BUFFER_OVERFLOW);

    for(int i=0; i < nvals; ++i) {
        val = StrStrMap_get(map, key[i], &status);
        assert(status==0);
        assert(0 == strcmp(val, key[i]));
    }

    // A key not in map
    val = StrStrMap_get(map, "not in map", &status);
    assert(status == 1);
    assert(0 == StrStrMap_hasKey(map, "not in map"));

    StrStrMap_free(map);

    unitTstResult("StrStrMap", "OK");

    return 0;
}
