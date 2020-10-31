#ifndef PtrU32MAP_H
#  define PtrU32MAP_H
#  include "typedefs.h"

#define MAPTYPE PtrU32Map
#define KEYTYPE void *
#define VALTYPE uint32_t
#define HASH uint64Hash
#define KEYCAST uint64_t

#define CMP(a,b) ((a)<(b) ? -1 : (a)>(b) ? 1 : 0)

#include "hashmap.h"

#endif
