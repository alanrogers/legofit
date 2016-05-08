#ifndef PATPROB_INCLUDED
#  define PATPROB_INCLUDED

#include "typedefs.h"

unsigned patprob(unsigned maxpat,
                 tipId_t pat[maxpat],
                 double prob[maxpat],
				 LblNdx *lblndx,
                 int nTasks,
				 long reps[nTasks],
                 int pointNdx,
                 const char *fname,
                 Bounds bnd);

#endif
