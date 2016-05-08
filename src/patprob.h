#ifndef PATPROB_INCLUDED
#  define PATPROB_INCLUDED

#include "typedefs.h"

unsigned patprob(unsigned maxpat,
                 tipId_t pat[maxpat],
                 double prob[maxpat],
                 int nTasks,
                 unsigned long nreps,
                 int pointNdx,
                 const char *fname,
                 Bounds bnd);

#endif
