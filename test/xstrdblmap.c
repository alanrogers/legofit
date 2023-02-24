/**
 * @file xstrdblmap.c
 * @author Alan R. Rogers
 * @brief Test strdblmap.c.
 * @copyright Copyright (c) 2023, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "strdblmap.h"
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
            fprintf(stderr, "usage: xstrdblmap [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xstrdblmap [-v]\n");
        exit(EXIT_FAILURE);
    }


    // Allocate small hash map so that it will be resized.
    StrDblMap *map = StrDblMap_new(8);

    assert(0 == StrDblMap_hasKey(map, "asdf"));

    const int nvals = 25;
    char *key[nvals];
    double value[nvals], val;
    int status;

    for(int i=0; i < nvals; ++i) {
        key[i] = malloc(2);
        key[i][0] = 'a' + i;
        key[i][1] = '\0';
        value[i] = (double) i;
        status = StrDblMap_insert(map, key[i], value[i]);
        assert(status == 0);
    }

    for(int i=0; i < nvals; ++i) {
        assert(1 == StrDblMap_hasKey(map, key[i]));
    }

    const char *longkey = "asdfasfadfasf asdfasfda fasdf adf asdf asdf ";
    status = StrDblMap_insert(map, longkey, 1.0);
    assert(status == 0);
    assert(1 == StrDblMap_hasKey(map, longkey));
    
    unsigned long size = StrDblMap_size(map);
    assert(nvals+1 == size);

    char * kk[size];
    status = StrDblMap_keys(map, size, kk);
    assert(status == 0);

    if(verbose) {
        for(int i=0; i<size; ++i) {
            val = StrDblMap_get(map, kk[i], &status);
            assert(status==0);
            printf("%s => %lf\n", kk[i], val);
        }                 
    }
    for(int i=0; i < size; ++i)
        free(kk[i]);

    char * kkk[size/2];
    status = StrDblMap_keys(map, size/2, kkk);
    assert(status == BUFFER_OVERFLOW);

    for(int i=0; i < nvals; ++i) {
        val = StrDblMap_get(map, key[i], &status);
        assert(status==0);
        assert(val == value[i]);
    }

    // A key not in map
    val = StrDblMap_get(map, "not in map", &status);
    assert(status == 1);
    assert(0 == StrDblMap_hasKey(map, "not in map"));

    StrDblMap_free(map);

    unitTstResult("StrDblMap", "OK");

    return 0;
}
