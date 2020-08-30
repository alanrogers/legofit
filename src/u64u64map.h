#ifndef U64U64MAP_H
#  define U64U64MAP_H
#  include "typedefs.h"
#  include <stdint.h>

#define MAPTYPE U64U64Map
#define KEYTYPE uint64_t
#define KEYCAST uint64_t
#define VALTYPE uint64_t
#define HASH uint64Hash
#define CMP(a,b) ((a)<(b) ? -1 : (a)>(b) ? 1 : 0)

#include "hashmap.h"

#endif
