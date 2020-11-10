#ifndef ARR_SEGMENT_H
#define ARR_SEGMENT_H

#include "typedefs.h"
#include <stdio.h>

int      Segment_addChild(void * vparent, void * vchild);
int      Segment_coalesce(Segment *self, int dosing, BranchTab *branchtab,
                          long unsigned *event_counter);
Segment *Segment_dup(Segment *old_root, PtrPtrMap *ppm);
int      Segment_equals(Segment *a, Segment *b);
int      Segment_feasible(const Segment * self, Bounds bnd, int verbose);
void     Segment_free(Segment *self);
void    *Segment_new(int twoN_i, int start_i, ParStore *ps);
void     Segment_newSample(Segment * self, unsigned ndx);
int      Segment_mix(void * vchild, int mix_i, void * vintrogressor,
                     void * vnative, ParStore *ps);
void     Segment_print(FILE * fp, void * self, int indent);
void     Segment_prune(Segment *self);
void    *Segment_root(void * vself);
void     Segment_sanityCheck(Segment * self, const char *file, int lineno);
void     Segment_unvisit(Segment *self);
void     Segment_update(Segment *self, ParStore *ps);
void     Segment_clear(Segment * self);
int      Segment_isClear(const Segment * self);

#endif
