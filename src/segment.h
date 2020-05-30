#ifndef ARR_SEGMENT_H
#define ARR_SEGMENT_H

#include "typedefs.h"
#include <stdio.h>

// One segment of a population network. This version works
// with MCTree.
struct Segment {
    int             nparents, nchildren, nsamples;
    double         *twoN;        // ptr to current pop size
    double         *start, *end; // duration of this PopNode
    double         *mix;         // ptr to frac of pop derived from parent[1]
    struct Segment *parent[2];
    struct Segment *child[2];

    tipId_t         sampId[MAXSAMP];

    int max;       // max number of lineages in segment

    // p[0][i] is prob there are i+1 lineages at recent end of segment
    // p[1][i] is analogous prob for ancient end of interval.
    double p[2][MAXSAMP];

    // Arrays of pointers to linked lists of IdSet objects. Dimension
    // is max X 2.  ids[0] refers to the recent end of the segment and
    // ids[1] to the ancient end. ids[0][i] is the list for the case
    // in which there are i+1 lineages at the recent end of the segment.
    IdSet *ids[2][MAXSAMP];
};

void    *Segment_new(double *twoN, double *start, NodeStore *ns);
int      Segment_coalesce(Segment *self, int maxsamp, int dosing,
                          BranchTab *branchtab, double v);


int      Segment_addChild(void * vparent, void * vchild);
int      Segment_mix(void * vchild, double *mPtr, void * vintrogressor, 
                      void * vnative);
void    *Segment_root(void * vself);
void     Segment_print(FILE * fp, void * self, int indent);

#endif
