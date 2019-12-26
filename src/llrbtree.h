#ifndef LEGOFIT_LLRBTREE_H
#define LEGOFIT_LLRBTREE_H

#include "typedefs.h"
#include <stdio.h>

typedef struct LlrbTree LlrbTree;

LlrbTree *LlrbTree_insert(LlrbTree *root, const void *key, const void *val,
                          void *(*dupKey)(const void *),
                          void *(*dupVal)(const void *),
                          int (*cmp)(const void *, const void *));
LlrbTree *LlrbTree_add(LlrbTree *root, const void *key, const void *val,
                       void *(*dupKey)(const void *),
                       void *(*dupVal)(const void *),
                       int (*cmp)(const void *, const void *),
                       void (*addVal)(void *, const void *));
void *LlrbTree_search(LlrbTree *root, const void *key,
                      int (*cmp)(const void *, const void *));
void Llbtree_free(LlrbTree *h);
void LlrbTree_free(LlrbTree *h, void (*freeKey)(void *),
                   void (*freeVal)(void *));
void LlrbTree_print(LlrbTree *h, FILE *fp, int indent,
                    int (*prKey)(const void *, FILE *fp),
                    int (*prVal)(const void *, FILE *fp));
#endif
