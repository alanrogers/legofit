#ifndef HAVE_PARSE
#define HAVE_PARSE

#include "typedefs.h"

int         countSegments(FILE * fp);
PtrPair     mktree(FILE * fp, SampNdx *sndx, LblNdx *lndx, Bounds *bnd,
                   NodeStore *ns);

#endif
