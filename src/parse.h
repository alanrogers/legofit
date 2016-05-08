#ifndef HAVE_PARSE
#define HAVE_PARSE

#include "typedefs.h"
PopNode    *mktree(FILE * fp, HashTab * poptbl, SampNdx *sndx, 
				   LblNdx *lndx, ParStore *parstore, Bounds *bnd);

#endif
