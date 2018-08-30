#ifdef CLASS
#  error "CLASS already defined"
#endif
#ifdef METHOD
#  error "METHOD already definedd"
#endif
#ifdef CAT
#  error "CAT already definedd"
#endif
#ifdef CAT1
#  error "CAT1 already definedd"
#endif

#define CLASS StrParMap
#define KEYTYPE const char *

#define CAT2(A,B) A ## _ ## B
#define CAT(A,B) CAT2(A, B)
#define METHOD(B) CAT(CLASS,B)

#include "parmap.h"

#undef CLASS
#undef METHOD
#undef CAT
#undef CAT1
