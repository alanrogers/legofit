#ifndef HAVE_PARSE
#define HAVE_PARSE

#include "typedefs.h"
int         countSegments(FILE * fp);
PopNode    *mktree(FILE * fp, SampNdx *sndx, LblNdx *lndx, ParStore *parstore,
                   Bounds *bnd, NodeStore *ns);

#endif
