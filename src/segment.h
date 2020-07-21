#ifndef ARR_SEGMENT_H
#define ARR_SEGMENT_H

#include "typedefs.h"
#include <stdio.h>

void    *Segment_new(int twoN_i, int start_i, ParStore *ps);
double   Segment_twoN(Segment *self);
double   Segment_start(Segment *self);
double   Segment_end(Segment *self);
double   Segment_mix(Segment *self);
int      Segment_coalesce(Segment *self, int maxsamp, int dosing,
                          BranchTab *branchtab);
int      Segment_addChild(void * vparent, void * vchild);
int      Segment_mix(void * vchild, double *mPtr, void * vintrogressor, 
                      void * vnative);
void     Segment_newSample(Segment * self, unsigned ndx);
void    *Segment_root(void * vself);
void     Segment_print(FILE * fp, void * self, int indent);
void     Segment_sanityCheck(Segment * self, const char *file, int lineno);
void     Segment_allocArrays(Segment *self, Stirling2 *stirling2);

#endif
