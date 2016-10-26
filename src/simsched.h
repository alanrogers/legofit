#ifndef ARR_SIMSCHED
#  define ARR_SIMSCHED

typedef struct SimSched SimSched;

SimSched   *SimSched_new(int nOptItr, int nSimReps);
void        SimSched_free(SimSched *self);
SimSched   *SimSched_append(SimSched * self, int nOptItr, int nSimReps);
void        SimSched_free(SimSched * self);
int         SimSched_getSimReps(SimSched * self);
int         SimSched_getOptItr(SimSched * self);
int         SimSched_next(SimSched * self);

#endif
