/**
 * @file llrbtree.c
 * @brief A left-leaning red-black tree.
 * 
 * These functions are generic, storing keys and values as void pointers.
 * They require as arguments function pointers that are able to interpret
 * the keys and values in terms of some concrete types.
 *
 * Translated from Java code at
 * http://www.cs.princeton.edu/~rs/talks/LLRB/LLRB.pdf 
 */
#include "llrbtree.h"
#include "strparmap.h"
#include "misc.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#define RED 1
#define BLACK 0
#define ISRED(node) ((node)!=NULL && (node)->color==RED)

struct LlrbTree {
    void *key;
    void *val;
    int color;
    LlrbTree *left, *right;
};

LlrbTree *LlrbTree_new(const void *key, const void *val,
                     void *(*dupKey)(const void *),
                     void *(*dupVal)(const void *));
LlrbTree *LlrbTree_insert_r(LlrbTree *h, const void *key, const void *val,
                          void *(*dupKey)(const void *),
                          void *(*dupVal)(const void *),
                          int (*cmp)(const void *, const void *));
LlrbTree *LlrbTree_add_r(LlrbTree *h, const void *key, const void *val,
                       void *(*dupKey)(const void *),
                       void *(*dupVal)(const void *),
                       int (*cmp)(const void *, const void *),
                       void (*addVal)(void *, const void *));
LlrbTree *LlrbTree_rotateLeft(LlrbTree *h);
LlrbTree *LlrbTree_rotateRight(LlrbTree *h);
void     LlrbTree_flipColors(LlrbTree *h);

void LlrbTree_free(LlrbTree *h, void (*freeKey)(void *),
                  void (*freeVal)(void *)) {
    if(h == NULL)
        return;
    LlrbTree_free(h->left, freeKey, freeVal);
    LlrbTree_free(h->right, freeKey, freeVal);
    (*freeKey)(h->key);
    (*freeVal)(h->val);
    free(h);
}

LlrbTree *LlrbTree_new(const void *key, const void *val,
                     void *(*dupKey)(const void *),
                     void *(*dupVal)(const void *)) {
    LlrbTree *self = malloc(sizeof(LlrbTree));
    CHECKMEM(self);
    self->key = (*dupKey)(key);
    self->val = (*dupVal)(val);
    self->color = RED;
    self->left = self->right = NULL;
    return self;
}

void *LlrbTree_search(LlrbTree *root, const void *key,
                     int (*cmp)(const void *, const void *)) {
    LlrbTree *x = root;
    while(x != NULL) {
        fprintf(stderr,"%s:%d\n",__FILE__,__LINE__);
        int diff = (*cmp)(key, x->val);
        if(diff == 0)
            return x->val;
        else if(diff < 0)
            x = x->left;
        else
            x = x->right;
    }
    return NULL;
}

LlrbTree *LlrbTree_insert(LlrbTree *root, const void *key, const void *val,
                        void *(*dupKey)(const void *),
                        void *(*dupVal)(const void *),
                        int (*cmp)(const void *, const void *)) {
    root = LlrbTree_insert_r(root, key, val, dupKey, dupVal,
                            cmp);
    root->color = BLACK;
    return root;
}

LlrbTree *LlrbTree_insert_r(LlrbTree *h, const void *key, const void *val,
                          void *(*dupKey)(const void *),
                          void *(*dupVal)(const void *),
                          int (*cmp)(const void *, const void *)) {
    if(h == NULL)
        return LlrbTree_new(key, val, dupKey, dupVal);
    if(ISRED(h->left) && ISRED(h->right))
        LlrbTree_flipColors(h);
    fprintf(stderr,"%s:%d key=%p h->key=%p\n",
            __FILE__,__LINE__, key, h->key);
    int diff = (*cmp)(key, h->key);
    fprintf(stderr,"%s:%d\n",__FILE__,__LINE__);
    if(diff == 0) {
        fprintf(stderr,"%s:%d: duplicate key\n", __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    if( diff < 0 )
        h->left = LlrbTree_insert_r(h->left, key, val, dupKey,
                                   dupVal, cmp);
    else {
        assert(diff > 0);
        h->right = LlrbTree_insert_r(h->right, key, val, dupKey,
                                    dupVal, cmp);
    }
    if(ISRED(h->right) && !ISRED(h->left))
        h = LlrbTree_rotateLeft(h);
    if(h->left && ISRED(h->left) && ISRED(h->left->left))
        h = LlrbTree_rotateRight(h);
    return h;
}

LlrbTree *LlrbTree_add(LlrbTree *root, const void *key, const void *val,
                     void *(*dupKey)(const void *),
                     void *(*dupVal)(const void *),
                     int (*cmp)(const void *, const void *),
                     void (*addVal)(void *, const void *)) {
    root = LlrbTree_add_r(root, key, val, dupKey, dupVal, cmp, addVal);
    root->color = BLACK;
    return root;
}

LlrbTree *LlrbTree_add_r(LlrbTree *h, const void *key, const void *val,
                       void *(*dupKey)(const void *),
                       void *(*dupVal)(const void *),
                       int (*cmp)(const void *, const void *),
                       void (*addVal)(void *, const void *)) {
    if(h == NULL)
        return LlrbTree_new(key, val, dupKey, dupVal);
    if(ISRED(h->left) && ISRED(h->right))
        LlrbTree_flipColors(h);
    fprintf(stderr,"%s:%d\n",__FILE__,__LINE__);
    int diff = (*cmp)(key, h->key);
    if(diff == 0)
        (*addVal)(h->val, val);
    else if( diff < 0 )
        h->left = LlrbTree_insert_r(h->left, key, val, dupKey,
                                   dupVal, cmp);
    else {
        assert(diff > 0);
        h->right = LlrbTree_insert_r(h->right, key, val, dupKey,
                                    dupVal, cmp);
    }
    if(ISRED(h->right) && !ISRED(h->left))
        h = LlrbTree_rotateLeft(h);
    if(h->left && ISRED(h->left) && ISRED(h->left->left))
        h = LlrbTree_rotateRight(h);
    return h;
}

LlrbTree *LlrbTree_rotateLeft(LlrbTree *h) {
    LlrbTree *x = h->right;
    h->right = x->left;
    x->left = h;
    x->color = h->color;
    h->color = RED;
    return x;
}

LlrbTree *LlrbTree_rotateRight(LlrbTree *h) {
    LlrbTree *x = h->left;
    h->left = x->right;
    x->right = h;
    x->color = h->color;
    h->color = RED;
    return x;
}

void LlrbTree_flipColors(LlrbTree *h) {
    h->color = !h->color;
    h->left->color = !h->left->color;
    h->right->color = !h->right->color;
}

void LlrbTree_print(LlrbTree *h, FILE *fp, int indent,
                    int (*prKey)(const void *, FILE *fp),
                    int (*prVal)(const void *, FILE *fp)) {
    if(h == NULL)
        return;
    LlrbTree_print(h->left, fp, indent+1, prKey, prVal);
    int i;
    for(i=0; i<indent; ++i)
        fprintf(fp,"%d", i);
    putc('[', fp);
    (*prKey)(h->key, fp);
    putc(',', fp);
    (*prVal)(h->val, fp);
    putc(',', fp);
    fprintf(fp, "%s]\n", h->color==RED ? "red" : "black");
    LlrbTree_print(h->right, fp, indent+1, prKey, prVal);
}
