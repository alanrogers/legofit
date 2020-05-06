/**
 * @file segment.c
 * @author Alan R. Rogers
 * @brief A single segment of a population tree.
 *
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "segment.h"
#include "idset.h"
#include "partprob.h"
#include "intpart.h"
#include "comb.h"
#include "binary.h"
#include "matcoal.h"
#include "misc.h"
#include "branchtab.h"
#include <stdlib.h>

typedef struct CombDat CombDat;

struct CombDat {
    double contribution; // Pr[site pattern]*E[len of interval]
    IdSet *ids;
    BranchTab *branchtab;
    int dosing;    // do singleton site patterns if nonzero
};

// Probabilistic variables for one segment of a population tree, i.e.,
// for ont PopNode.
struct Segment {
    int allocated; // allocated dimension of arrays
    int max;     // max number of lineages in segment

    // p[0][i] is prob there are i+1 lineages at recent end of segment
    // p[1][i] is analogous prob for ancient end of interval.
    double *p[2];

    // Arrays of pointers to linked lists of IdSet objects. Dimension
    // is max X 2.  ids[0] refers to the recent end of the segment and
    // ids[1] to the ancient end. ids[0][i] is the list for the case
    // in which there are i+1 lineages at the recent end of the segment.
    IdSet **ids[2];
};

int    visitComb(int d, int ndx[d], void *data);

Segment *Segment_new(void) {
    Segment *self = malloc(sizeof(Segment));
    CHECKMEM(self);

    self->allocated = 8;
    self->max = 0;

    self->p[0] = malloc(self->allocated * sizeof(self->p[0][0]));
    CHECKMEM(self->p[0]);

    self->p[1] = malloc(self->allocated * sizeof(self->p[1][0]));
    CHECKMEM(self->p[1]);

    self->ids[0] = self->ids[1] = NULL;

    return self;
}

// Add IdSet to Segment
void Segment_add(Segment *self, IdSet *idset) {
    int nids = idset->nids;

    // reallocate arrays if necessary
    if(nids > self->allocated) {
        if(2*self->allocated >= nids)
            self->allocated *= 2;
        else
            self->allocated = nids;

        self->p[0] = realloc(self->p[0],
                             self->allocated * sizeof(self->p[0][0]));
        CHECKMEM(self->p[0]);
        self->p[1] = realloc(self->p[1],
                             self->allocated * sizeof(self->p[1][0]));
        CHECKMEM(self->p[1]);

        self->ids[0] = realloc(self->ids[0],
                               self->allocated * sizeof(self->ids[0][0]));
        CHECKMEM(self->ids[0]);
        self->ids[1] = realloc(self->ids[1],
                               self->allocated * sizeof(self->ids[1][0]));
        CHECKMEM(self->ids[1]);
    }

    // add IdSet to segment
    self->ids[0][nids-1] = IdSet_add(self->ids[0][nids-1], nids,
                                     idset->tid, 0.0);
}

int Segment_coalesce(Segment *self, int maxsamp, int dosing,
                     BranchTab *branchtab, double v) {
    assert(self->max <= maxsamp);

    double pr[maxsamp], elen[maxsamp], sum;
    int n, i, status=0;
    const int finite = isfinite(v); // is segment finite?
    int kstart = (finite ? 1 : 2);  // skip k=1 if infinite

    // If there is only one line of descent, no coalescent events are
    // possible, so p[1][0] is at least as large as p[0][0].
    self->p[1][0] = self->p[0][0];

    // Calculate probabilities and expected values, p[1] and elen.
    for(n=2; n <= self->max; ++n) {

        // Skip improbable states.
        if(self->p[0][n-1] == 0.0)
            continue;
        
        // Calculate prob of 2,3,...,n lines of descent.
        // On return, pr[1] = prob[2], pr[n-1] = prob[n],
        // pr[0] not set.
        MatCoal_project(n-1, pr+1, v);

        // Calculate probability of 1 line of descent. I'm doing the
        // sum in reverse order on the assumption that higher indices
        // will often have smaller probabilities.
        // pr[i] is prob of i+1 lineages; pr[0] not set
        sum=0.0;
        for(i=n; i > 1; --i)
            sum += pr[i-1];
        self->p[1][0] += self->p[0][n-1] * (1.0-sum);

        // Add probs of 2..n lines of descent
        // pr[i] is prob of i+2 lineages
        // self->p[1][i] is prob if i+1 lineages
        for(i=2; i <= n; ++i)
            self->p[1][i-1] += self->p[0][n-1] * pr[i-1];

        // Calculate expected length, within the segment, of
        // coalescent intervals with 2,3,...,n lineages.  elen[i] is
        // expected length of interval with i-1 lineages.
        // elen[0] not set.
        MatCoal_ciLen(n-1, elen+1, v);

        // Calculate expected length, elen[0], of interval with one
        // lineage.  elen[i] refers to interval with i+1 lineages.
        sum = 0.0;
        for(i=n; i >= 2; --i)
            sum += elen[i-1];
        elen[0] = v - sum; // elen[0] is infinite if v is.

        // Multiply elen by probability that segment had n lineages on
        // recent end of segment.
        for(i=1; i <= n; ++i)
            elen[i-1] *= self->p[0][n-1];
    }

    CombDat cd = {.contribution = 0.0,
                  .ids = NULL,
                  .branchtab = branchtab,
                  .dosing = dosing
    };

    // loop over number of descendants in this segment
    for(n=1; n <= self->max; ++n) {
        // Skip improbable states.
        if(self->p[0][n-1] == 0.0)
            continue;

        cd.ids = ids[0][n-1];

        // loop over intervals: k is the number of ancestors
        // w/i interval.
        for(int k=kstart; k <= n; ++k) {
            
            // portion of log Qdk that doesn't involve d
            long double lnconst = logl(k) - lbinom(n-1, k-1);

            // Within each interval, there can be ancestors
            // with 1 descendant, 2, 3, ..., n-k+2.
            for(int d=1; d <= n-k+1; ++d) {
                long double lnprob = lnconst
                    + lbinom(n-d-1, k-1) - lbinom(n,d);

                // probability of site pattern
                cd.contribution = (double) expl(lnprob);

                // times expected length of interval
                cd.contribution *= elen[k-1];
                
                status = traverseComb(n, d, visitComb, &cd);
                if(status)
                    return status;
            }
        }
    }
    return status;
}

/// Visit a combination
int visitComb(int d, int ndx[d], void *data) {
    assert(d>0);
    CombDat *dat = (CombDat *) data;

    for(IdSet ids = dat->ids; ids != NULL; ids = ids->next) {
        
        // tid is the union of descendants in current set.
        tipId_t tid = 0;
        for(int i=0; i < d; ++i)
            tid |= ids->tid[ndx[i]];

        // Skip singletons unless data->dosing is nonzero
        if(!dat->dosing && isPow2(tid[i]))
            continue;
      
        // Increment BranchTab entry for current tid value.
        BranchTab_add(dat->branchtab, tid,
                      ids->p * dat->contribution);
    }
    return 0;
}
