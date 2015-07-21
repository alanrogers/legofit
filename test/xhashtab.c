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
#include "gptree.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

void prval(void *vx, void *vfp);

// For use by HashTab_map. It will print every node in the table.
void prval(void *vx, void *vfp) {
    PopNode *node = (PopNode *) vx;
    FILE *fp = (FILE *) vfp;
    PopNode_printShallow(node, fp);
}

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

    printf("Calling HashTab_print\n");
    HashTab_print(ht);

    printf("Calling HashTab_map to print table\n");
    HashTab_map(ht, prval, stdout);

    printf("Number of elements in HashTab: %lu\n",
           HashTab_size(ht));

    HashTabSeq *hts = HashTabSeq_new(ht);
    CHECKMEM(hts);
    e = HashTabSeq_next(hts);
    while(e != NULL) {
        El_print(e);
        putchar('\n');
        e = HashTabSeq_next(hts);
    }

    HashTabSeq_free(hts);
    HashTab_free(ht);
    free(key1);
    free(val1);
    free(key2);
    free(val2);
    free(key3);
    free(val3);

    return 0;
}
