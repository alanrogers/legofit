#ifndef ARR_GPTREE_H
#  define ARR_GPTREE_H

#  include "typedefs.h"
#  include <gsl/gsl_rng.h>


GPTree     *GPTree_new(const char *fname, Bounds bnd);
void        GPTree_free(GPTree *self);

NodeStore  *NodeStore_new(int len, PopNode *v);
void        NodeStore_free(NodeStore *self);
PopNode    *NodeStore_alloc(NodeStore *self);

#endif
