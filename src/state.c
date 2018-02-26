#include "state.h"
#include <assert.h>
#include <stdlib.h>

typedef struct State State;

struct State {
    int npts, npar; // numbers of points and parameters
    double **s;     // s[i][j]=value of j'th param at i'th point
};

// Allocate a new State object. There is no point in optimizing
// this, because State is used only at the beginning and end
// of Legofit.
State *State_new(int npts, int npar) {
    int i;
    State *self = malloc(sizeof(State));
    CHECKMEM(self);
    self->npts = npts;
    self->npar = npar;
    self->s = malloc(npts * sizeof(self->s[0]));
    CHECKMEM(self->s);
    for(i=0; i < npar; ++i) {
        self->s[i] = malloc(npar * sizeof(self->s[0][0]));
        CHECKMEM(self->s[i]);
    }
    return self;
}

void State_free(State *self) {
    assert(self);
    int i;
    for(i=0; i < self->npts; ++i)
        free(self->s[i]);
    free(self->s);
    free(self);
}
