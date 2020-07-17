#ifndef PtrPtrMAP_H
#  define PtrPtrMAP_H
#  include "typedefs.h"

#define MAPTYPE PtrPtrMap
#define KEYTYPE void *
#define VALTYPE void *
#define HASH uint64Hash
#define KEYCAST uint64_t

#define HASH_DIM 256ul
#define CMP(a,b) ((a)<(b) ? -1 : (a)>(b) ? 1 : 0)

#include "hashmap.h"

#endif
