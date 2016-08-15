#ifndef HAVE_TYPEDEFS
#  define HAVE_TYPEDEFS

#include <stdint.h>
#define FILENAMESIZE 200

typedef struct Boot Boot;
typedef struct BootChr BootChr;
typedef struct Bounds Bounds;
typedef struct BranchTab BranchTab;
typedef struct El El;
typedef struct Gene Gene;
typedef struct GPTree GPTree;
typedef struct HashTab HashTab;
typedef struct HashTabSeq HashTabSeq;
typedef struct LblNdx LblNdx;
typedef struct NodeStore NodeStore;
typedef struct ParKeyVal ParKeyVal;
typedef struct ParStore ParStore;
typedef struct PopNode PopNode;
typedef struct SampNdx SampNdx;
typedef struct StrInt StrInt;
typedef struct Tokenizer Tokenizer;
typedef struct DAFReader DAFReader;

#if 1
#  define TIPID_SIZE 32
typedef uint32_t tipId_t;
#else
#  define TIPID_SIZE 64
typedef uint64_t tipId_t;
#endif

#endif
