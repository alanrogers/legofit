#ifndef U64I64MAP_H
#  define U64I64MAP_H
#  include "typedefs.h"
#  include <stdint.h>

#define MAPTYPE U64I64Map
#define KEYTYPE uint64_t
#define KEYCAST uint64_t
#define VALTYPE int64_t
#define HASH uint64Hash
#define HASH_DIM 1024ul
#define CMP(a,b) ((a)<(b) ? -1 : (a)>(b) ? 1 : 0)

#include "hashmap.h"

#endif
