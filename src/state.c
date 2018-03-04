#include "state.h"
#include "error.h"
#include "misc.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

struct State {
    int npts, npar; // numbers of points and parameters
    double **s;     // s[i][j]=value of j'th param at i'th point
};

int State_npoints(State *self) {
    return self->npts;
}

int State_nparameters(State *self) {
    return self->npar;
}

// Set state vector with index "ndx" equal to vector x.
int State_setVector(State *self, int ndx, int dim, double x[dim]) {
    if(dim != self->npar)
        return EINVAL;
    if(ndx >= self->npts)
        return EINVAL;
    memcpy(self->s[ndx], x, dim * sizeof(x[0]));
    return 0;
}

// Copy state vector with index "ndx" into vector x.
int State_getVector(State *self, int ndx, int dim, double x[dim]) {
    if(dim != self->npar)
        return EINVAL;
    if(ndx >= self->npts)
        return EINVAL;
    memcpy(x, self->s[ndx], self->npar * sizeof(self->s[0]));
    return 0;
}

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
    for(i=0; i < npts; ++i) {
        self->s[i] = malloc(npar * sizeof(double));
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

// Construct a new State object by reading a file
State *State_read(FILE *fp) {
    int i, j, npts, npar, status;
    status = fscanf(fp, "%d %d", &npts, &npar);
    if(status != 2) {
        fprintf(stderr,"%s:%d: Can't read dimensions in state file\n",
                __FILE__,__LINE__);
        return NULL;
    }
    State *self = State_new(npts, npar);
    CHECKMEM(self);

    for(i=0; i < npts; ++i) {
        for(j=0; j < npar; ++j) {
            status = fscanf(fp, "%lf", self->s[i]+j);
            if(status != 1) {
                fprintf(stderr,"%s:%d:"
                        " Can't read value (%d,%d) in state file\n",
                        __FILE__,__LINE__,i,j);
                State_free(self);
                return NULL;
            }
        }
    }
    return self;
}

// Print State object to a file
int State_print(State *self, FILE *fp) {
    int i, j, status;

    status = fprintf(fp, "%d %d\n", self->npts, self->npar);
    if(status==0) {
        fprintf(stderr,"%s:%d: can't write to file\n",
                __FILE__,__LINE__);
        return EIO;
    }
    for(i=0; i < self->npts; ++i) {
        for(j=0; j < self->npar; ++j) {
            status = fprintf(fp, " %0.18lf", self->s[i][j]);
            if(status == 0) {
                fprintf(stderr,"%s:%d: can't write to file\n",
                        __FILE__,__LINE__);
                return EIO;
            }
        }
        putc('\n', fp);
    }
    return 0;
}
