// Translated from Java code at
// http://www.cs.princeton.edu/~rs/talks/LLRB/LLRB.pdf
#include "llrbtree.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#define RED 1
#define BLACK 0
#define ISRED(node) ((node)!=NULL && (node)->color==RED)

struct BSTNode {
    bstkey_t key;
    val_t value;
    int color;
    BSTNode *left, *right;
};

BSTNode *BSTNode_new(bstkey_t key, val_t value);
BSTNode *BSTNode_insert_r(BSTNode *h, bstkey_t key, val_t value);
BSTNode *BSTNode_rotateLeft(BSTNode *h);
BSTNode *BSTNode_rotateRight(BSTNode *h);
void BSTNode_flipColors(BSTNode *h);

void BSTNode_free(BSTNode *h) {
    if(h == NULL)
        return;
    BSTNode_free(h->left);
    BSTNode_free(h->right);
    free(h);
}

val_t BSTNode_value(BSTNode *h) {
    return h->value;
}

bstkey_t BSTNode_key(BSTNode *h) {
    return h->key;
}

BSTNode *BSTNode_new(bstkey_t key, val_t value) {
    BSTNode *self = malloc(sizeof(BSTNode));
    if(self==NULL) {
        fprintf(stderr,"%s:%d: can't allocate memory\n",
                __FILE__, __LINE__);
        exit(1);
    }
    self->key = key;
    self->value = value;
    self->color = RED;
    return self;
}

BSTNode *BSTNode_search(BSTNode *root, bstkey_t key) {
    BSTNode *x = root;
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

BSTNode *BSTNode_insert(BSTNode *root, bstkey_t key, val_t value) {
    root = BSTNode_insert_r(root, key, value);
    root->color = BLACK;
    return root;
}

BSTNode *BSTNode_insert_r(BSTNode *h, bstkey_t key, val_t value) {
    if(h == NULL)
        return BSTNode_new(key, value);
    if(ISRED(h->left) && ISRED(h->right))
        BSTNode_flipColors(h);
    if(key == h->key)
        h->value = value;
    else if( key < h->key )
        h->left = BSTNode_insert_r(h->left, key, value);
    else {
        assert(key > h->key);
        h->right = BSTNode_insert_r(h->right, key, value);
    }
    if(ISRED(h->right) && !ISRED(h->left))
        h = BSTNode_rotateLeft(h);
    if(ISRED(h->left) && ISRED(h->left->left))
        h = BSTNode_rotateRight(h);
    return h;
}

BSTNode *BSTNode_rotateLeft(BSTNode *h) {
    BSTNode *x = h->right;
    h->right = x->left;
    x->left = h;
    x->color = h->color;
    h->color = RED;
    return x;
}

BSTNode *BSTNode_rotateRight(BSTNode *h) {
    BSTNode *x = h->left;
    h->left = x->right;
    x->right = h;
    x->color = h->color;
    h->color = RED;
    return x;
}

void BSTNode_flipColors(BSTNode *h) {
    h->color = !h->color;
    h->left->color = !h->left->color;
    h->right->color = !h->right->color;
}

#ifdef TEST

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

#include <string.h>

void BSTNode_print(BSTNode *h, FILE *fp, int indent);

void BSTNode_print(BSTNode *h, FILE *fp, int indent) {
    if(h == NULL)
        return;
    BSTNode_print(h->left, fp, indent+1);
    int i;
    for(i=0; i<indent; ++i)
        fprintf(fp,"%d", i);

    // This fprintf statement must be modified to reflect
    // the definitions of bstkey_t and val_t.
    fprintf(fp, "[%u, %lf, %s]\n",
            h->key,
            h->value,
            h->color==RED ? "red" : "black");
    BSTNode_print(h->right, fp, indent+1);
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

    BSTNode *root = NULL;
    int i;
    for(i=0; i < 100; ++i) {
        bstkey_t key = (bstkey_t) i;
        val_t value = (val_t) i;
        root = BSTNode_insert(root, key, value);
    }

    if(verbose)
        BSTNode_print(root, stdout, 0);

    // Search for nodes that exist. Each search should succeed.
    for(i=0; i < 100; ++i) {
        bstkey_t key = (bstkey_t) i;
        val_t value = (val_t) i;
        BSTNode *found = BSTNode_search(root, key);
        assert(found);
        assert(BSTNode_key(found) == key);
        assert(BSTNode_value(found) == value);
    }

    // Search for nodes that don't exist. Each search should fail.
    for(i=100; i < 110; ++i) {
        bstkey_t key = (bstkey_t) i;
        BSTNode *found = BSTNode_search(root, key);
        assert(found == NULL);
    }
    BSTNode_free(root);
    printf("%-26s %s\n", "BSTNode", "OK");
}
#endif
