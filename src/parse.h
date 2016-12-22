#ifndef HAVE_PARSE
#define HAVE_PARSE

#include "typedefs.h"
int         countSegments(FILE * fp);
PopNode    *mktree(FILE * fp, SampNdx *sndx, LblNdx *lndx, ParStore *parstore,
                   ExoPar *ep, Bounds *bnd, NodeStore *ns);

#endif
