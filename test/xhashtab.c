/**
 * @file xhashtab.c
 * @author Alan R. Rogers
 * @brief Test hashtab.c.
 * @copyright Copyright (c) 2015, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "hashtab.h"
#include "misc.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(void) {
    El *e;
    HashTab *ht = HashTab_new();
    CHECKMEM(ht);
    assert(0 == HashTab_size(ht));

    char *key1 = strdup("key1");
    char *val1 = strdup("val1");
    e = HashTab_get(ht, key1);
    assert(1 == HashTab_size(ht));
    assert(NULL == El_get(e));
    El_set(e, val1);

    char *key2 = strdup("key2");
    char *val2 = strdup("val2");
    e = HashTab_get(ht, key2);
    assert(2 == HashTab_size(ht));
    assert(NULL == El_get(e));
    El_set(e, val2);

    char *key3 = strdup("key3");
    char *val3 = strdup("val3");
    e = HashTab_get(ht, key3);
    assert(3 == HashTab_size(ht));
    assert(NULL == El_get(e));
    El_set(e, val3);

    assert(0 == strcmp(val1, (char*) El_get(HashTab_get(ht, key1))));
    assert(0 == strcmp(val2, (char*) El_get(HashTab_get(ht, key2))));
    assert(0 == strcmp(val3, (char*) El_get(HashTab_get(ht, key3))));

    HashTab_print(ht);

    printf("Number of elements in HashTab: %lu\n",
           HashTab_size(ht));

    HashTab_free(ht);
    free(key1);
    free(val1);
    free(key2);
    free(val2);
    free(key3);
    free(val3);

    return 0;
}
