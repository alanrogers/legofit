#ifndef PATPROB_INCLUDED
#  define PATPROB_INCLUDED

#include "typedefs.h"
unsigned patprob(unsigned maxpat,
                 tipId_t pat[maxpat],
                 double prob[maxpat],
                 GPTree *gptree,
                 LblNdx *lblndx,
                 int nTasks,
                 long reps[nTasks],
                 int pointNdx);
#endif