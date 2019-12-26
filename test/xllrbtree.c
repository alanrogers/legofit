/**
 * @file xllrbtree.c
 * @author Alan R. Rogers
 * @brief Test llrbtree, a left-leaning red-black tree.
 *
 * This unit test assumes that keys are character strings and values
 * are ints.
 *
 * @copyright Copyright (c) 2019, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "llrbtree.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

#include <string.h>

void *dupKey(const void * key);
void *dupVal(const void * val);
int cmp(const void *x, const void *y);
void addVal(void *x, const void *y);
void freeKey(void *key);
void freeVal(void *val);
int prKey(const void *key, FILE *fp);
int prVal(const void *val, FILE *fp);

// duplicate a string
void *dupKey(const void * key) {
    return (void *) strdup( (const char *) key);
}

// duplicate an int
void *dupVal(const void * val) {
    int *dup = malloc(sizeof(*dup));
    CHECKMEM(dup);
    memcpy(dup, val, sizeof(*dup));
    return dup;
}

// compare strings
int cmp(const void *x, const void *y) {
    const char *sx = (const char *) x;
    const char *sy = (const char *) y;
    return strcmp(sx, sy);
}

// add right argument to left argument, interpreting as int pointers
void addVal(void *x, const void *y) {
    int *ix = (int *) x;
    const int *iy = (const int *) y;

    *ix += *iy;
}

void freeKey(void *key) {
    free(key);
}

void freeVal(void *val) {
    free(val);
}

// print a string
int prKey(const void *key, FILE *fp) {
    return fprintf(fp, "%s", (const char *) key);
}

// print an int
int prVal(const void *val, FILE *fp) {
    return fprintf(fp, "%d", *((const int *) val));
}

int main(int argc, char **argv) {

    int         verbose = 0;

    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xllrbtree [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    const char *key[] = {"key0", "key1", "key2", "key3"};
    int val[] = {0, 1, 2, 5};

    fprintf(stderr,"%s:%d\n",__FILE__,__LINE__);

    LlrbTree *root = NULL;
    root = LlrbTree_insert(root, key[0], val+0, dupKey, dupVal, cmp);
    root = LlrbTree_insert(root, key[1], val+1, dupKey, dupVal, cmp);
    root = LlrbTree_insert(root, key[2], val+2, dupKey, dupVal, cmp);
    root = LlrbTree_insert(root, key[3], val+3, dupKey, dupVal, cmp);

    fprintf(stderr,"%s:%d\n",__FILE__,__LINE__);

    if(verbose)
        LlrbTree_print(root, stderr, 0, prKey, prVal);

    // Search for nodes that exist. Each search should succeed.
    int *vptr = NULL;

    fprintf(stderr,"%s:%d: searching %s\n",__FILE__,__LINE__, key[0]);
    vptr = LlrbTree_search(root, key[0], cmp);
    assert(vptr);
    assert(*vptr == val[0]);

    vptr = LlrbTree_search(root, key[1], cmp);
    assert(vptr);
    assert(*vptr == val[1]);

    vptr = LlrbTree_search(root, key[2], cmp);
    assert(vptr);
    assert(*vptr == val[2]);

    vptr = LlrbTree_search(root, key[3], cmp);
    assert(vptr);
    assert(*vptr == val[3]);

    // Search for nodes that don't exist. Each search should fail.
    vptr = LlrbTree_search(root, "notthere", cmp);
    assert(vptr == NULL);

    // Add to a value
    int i=2;
    root = LlrbTree_add(root, key[1], &i, dupKey, dupVal, cmp, addVal);
    assert(root);
    vptr = LlrbTree_search(root, key[1], cmp);
    assert(vptr);
    assert(*vptr == val[1]+2);

    LlrbTree_free(root, freeKey, freeVal);

    unitTstResult("LlrbTree", "OK");
}
