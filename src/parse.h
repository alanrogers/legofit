#ifndef HAVE_PARSE
#define HAVE_PARSE

#include "typedefs.h"

PopNode    *mktree(FILE * fp, HashTab * ht, SampNdx *sndx, ParStore *fixed,
				   Bounds *bnd);

#endif
