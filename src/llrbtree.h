#include "typedefs.h"
typedef struct BSTNode BSTNode;
typedef tipId_t bstkey_t;
typedef double val_t;

val_t BSTNode_value(BSTNode *h);
bstkey_t BSTNode_key(BSTNode *h);
BSTNode *BSTNode_search(BSTNode *root, bstkey_t key);
BSTNode *BSTNode_insert(BSTNode *root, bstkey_t key, val_t value);
void BSTNode_free(BSTNode *h);
