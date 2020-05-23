#ifndef ARR_NODESTORE_H
#  define ARR_NODESTORE_H

#  include "typedefs.h"

NodeStore  *NodeStore_new(unsigned len, size_t elsize, void * v);
void        NodeStore_free(NodeStore *self);
void       *NodeStore_alloc(NodeStore * self);

#endif
