#ifndef HAVE_PARSE
#define HAVE_PARSE

#include "typedefs.h"

PopNode    *mktree(FILE * fp, HashTab * poptbl, SampNdx *sndx, ParStore *fixed,
				   ParStore *var, Bounds *bnd);

#endif
