#ifndef ARR_GPTREE_H
#  define ARR_GPTREE_H

#  include "typedefs.h"

GPTree     *GPTree_new(const char *fname, Bounds bnd);
void        GPTree_free(GPTree *self);
GPTree     *GPTree_dup(GPTree *old);

#endif
