#ifndef ARR_SIMSCHED
#  define ARR_SIMSCHED

#  include "typedefs.h"
#  include <stdio.h>

SimSched   *SimSched_new(void);
SimSched   *SimSched_dup(const SimSched *self);
int         SimSched_nStages(const SimSched *self);
void        SimSched_free(SimSched *self);
void        SimSched_append(SimSched * self, long nOptItr, long nSimReps);
void        SimSched_free(SimSched * self);
long        SimSched_getOptItr(SimSched * self);
long        SimSched_getSimReps(SimSched * self);
int         SimSched_next(SimSched * self);
void        SimSched_print(const SimSched *self, FILE *fp);

#endif
