/**
 * @file addrparmap.c
 * @brief Map parameter addresses to parameter structures.
 * 
 * This class uses a left-leaning binary tree to build a map
 * that associates parameter addresses with structures that include
 * a pointer to the current value, the lower and upper bounds, and
 * the type of the parameter. The interface allows parameters to be
 * inserted and later looked up by address.
 *
 * Translated from Java code at
 * http://www.cs.princeton.edu/~rs/talks/LLRB/LLRB.pdf 
 */
#include "param.h"
#include "strparmap.h"
#include "misc.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#define RED 1
#define BLACK 0
#define ISRED(node) ((node)!=NULL && (node)->color==RED)

struct AddrParMap {
    Param *par; // locally owned
    int color;
    AddrParMap *left, *right;
};

AddrParMap *AddrParMap_new(Param *par);
AddrParMap *AddrParMap_insert_r(AddrParMap *h, Param *par);
AddrParMap *AddrParMap_rotateLeft(AddrParMap *h);
AddrParMap *AddrParMap_rotateRight(AddrParMap *h);
void AddrParMap_flipColors(AddrParMap *h);

void AddrParMap_free(AddrParMap *h) {
    if(h == NULL)
        return;
    AddrParMap_free(h->left);
    AddrParMap_free(h->right);
    Param_free(h->par);
    free(h);
}

AddrParMap *AddrParMap_new(Param *par) {
    AddrParMap *self = malloc(sizeof(AddrParMap));
    if(self==NULL) {
        fprintf(stderr,"%s:%d: can't allocate memory\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    self->key = key;
    self->value = value;
    self->color = RED;
    return self;
}

AddrParMap *AddrParMap_search(AddrParMap *root, bstkey_t key) {
    AddrParMap *x = root;
    while(x != NULL) {
        if(key == x->key)
            return x;
        else if(key < x->key)
            x = x->left;
        else
            x = x->right;
    }
    return NULL;
}

AddrParMap *AddrParMap_insert(AddrParMap *root, bstkey_t key, val_t value) {
    root = AddrParMap_insert_r(root, key, value);
    root->color = BLACK;
    return root;
}

AddrParMap *AddrParMap_insert_r(AddrParMap *h, bstkey_t key, val_t value) {
    if(h == NULL)
        return AddrParMap_new(key, value);
    if(ISRED(h->left) && ISRED(h->right))
        AddrParMap_flipColors(h);
    if(key == h->key)
        h->value = value;
    else if( key < h->key )
        h->left = AddrParMap_insert_r(h->left, key, value);
    else {
        assert(key > h->key);
        h->right = AddrParMap_insert_r(h->right, key, value);
    }
    if(ISRED(h->right) && !ISRED(h->left))
        h = AddrParMap_rotateLeft(h);
    if(ISRED(h->left) && ISRED(h->left->left))
        h = AddrParMap_rotateRight(h);
    return h;
}

AddrParMap *AddrParMap_rotateLeft(AddrParMap *h) {
    AddrParMap *x = h->right;
    h->right = x->left;
    x->left = h;
    x->color = h->color;
    h->color = RED;
    return x;
}

AddrParMap *AddrParMap_rotateRight(AddrParMap *h) {
    AddrParMap *x = h->left;
    h->left = x->right;
    x->right = h;
    x->color = h->color;
    h->color = RED;
    return x;
}

void AddrParMap_flipColors(AddrParMap *h) {
    h->color = !h->color;
    h->left->color = !h->left->color;
    h->right->color = !h->right->color;
}

#ifdef TEST

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

#include <string.h>

void AddrParMap_print(AddrParMap *h, FILE *fp, int indent);

void AddrParMap_print(AddrParMap *h, FILE *fp, int indent) {
    if(h == NULL)
        return;
    AddrParMap_print(h->left, fp, indent+1);
    int i;
    for(i=0; i<indent; ++i)
        fprintf(fp,"%d", i);

    // This fprintf statement must be modified to reflect
    // the definitions of bstkey_t and val_t.
    fprintf(fp, "[%u, %lf, %s]\n",
            h->key,
            h->value,
            h->color==RED ? "red" : "black");
    AddrParMap_print(h->right, fp, indent+1);
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

    AddrParMap *root = NULL;
    int i;
    for(i=0; i < 100; ++i) {
        bstkey_t key = (bstkey_t) i;
        val_t value = (val_t) i;
        root = AddrParMap_insert(root, key, value);
    }

    if(verbose)
        AddrParMap_print(root, stdout, 0);

    // Search for nodes that exist. Each search should succeed.
    for(i=0; i < 100; ++i) {
        bstkey_t key = (bstkey_t) i;
        val_t value = (val_t) i;
        AddrParMap *found = AddrParMap_search(root, key);
        assert(found);
        assert(AddrParMap_key(found) == key);
        assert(AddrParMap_value(found) == value);
    }

    // Search for nodes that don't exist. Each search should fail.
    for(i=100; i < 110; ++i) {
        bstkey_t key = (bstkey_t) i;
        AddrParMap *found = AddrParMap_search(root, key);
        assert(found == NULL);
    }
    AddrParMap_free(root);
    printf("%-26s %s\n", "AddrParMap", "OK");
}
#endif
