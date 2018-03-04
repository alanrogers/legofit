#include "state.h"
#include "error.h"
#include "misc.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

struct State {
    int npts, npar; // numbers of points and parameters
    double *cost;   // cost[i] is cost function at i'th point
    double **s;     // s[i][j]=value of j'th param at i'th point
};

int State_npoints(State *self) {
    return self->npts;
}

int State_nparameters(State *self) {
    return self->npar;
}

void State_setCost(State *self, int ndx, double cost) {
    if(ndx >= self->npts) {
        fprintf(stderr,"%s:%d: index out of bounds\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    self->cost[ndx] = cost;
}

double State_getCost(State *self, int ndx) {
    if(ndx >= self->npts) {
        fprintf(stderr,"%s:%d: index out of bounds\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    return self->cost[ndx];
}


// Set state vector with index "ndx" equal to vector x.
void State_setVector(State *self, int ndx, int dim, double x[dim]) {
    if(dim != self->npar || ndx >= self->npts) {
        fprintf(stderr,"%s:%d: index out of bounds\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    memcpy(self->s[ndx], x, dim * sizeof(x[0]));
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
    self->cost = malloc(npts * sizeof(self->cost[0]));
    CHECKMEM(self->cost);
    self->s = malloc(npts * sizeof(self->s[0]));
    CHECKMEM(self->s);
    for(i=0; i < npts; ++i) {
        self->cost[i] = strtod("NaN", NULL);
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
    free(self->cost);
    free(self->s);
    free(self);
}

// Construct a new State object by reading a file
State *State_read(FILE *fp) {
    int i, j, npts, npar, status;
    State *self = NULL;
    status = fscanf(fp, "%d %d", &npts, &npar);
    if(status != 2) {
        fprintf(stderr,"%s:%d: status=%d\n", __FILE__,__LINE__,status);
        goto fail;
    }
    self = State_new(npts, npar);
    CHECKMEM(self);

    for(i=0; i < npts; ++i) {
        status = fscanf(fp, "%lf", self->cost + i);
        if(status != 1) {
            fprintf(stderr,"%s:%d: status=%d\n", __FILE__,__LINE__,status);
            goto fail;
        }
        for(j=0; j < npar; ++j) {
            status = fscanf(fp, "%lf", self->s[i]+j);
            if(status != 1) {
                fprintf(stderr,"%s:%d: status=%d\n", __FILE__,__LINE__,status);
                goto fail;
            }
        }
    }
    return self;

 fail:
    fprintf(stderr,"%s:%d: Can't read state file\n",
            __FILE__,__LINE__);
    if(self)
        State_free(self);
    return NULL;
}

// Print State object to a file
int State_print(State *self, FILE *fp) {
    int i, j, status, imin=0;

    // find index of minimum cost
    for(i=1; i < self->npts; ++i)
        if(self->cost[i] < self->cost[imin])
            imin = i;

    status = fprintf(fp, "%d %d\n", self->npts, self->npar);
    if(status==0)
        goto fail;

    // print point with minimum cost first
    status = fprintf(fp, "%0.18lg", self->cost[imin]);
    if(status == 0)
        goto fail;
    for(j=0; j < self->npar; ++j) {
        status = fprintf(fp, " %0.18lg", self->s[imin][j]);
        if(status == 0)
            goto fail;
    }
    putc('\n', fp);

    // print remainint lines
    for(i=0; i < self->npts; ++i) {
        if(i == imin)
            continue;
        status = fprintf(fp, "%0.18lg", self->cost[i]);
        if(status == 0)
            goto fail;
        for(j=0; j < self->npar; ++j) {
            status = fprintf(fp, " %0.18lg", self->s[i][j]);
            if(status == 0)
                goto fail;
        }
        putc('\n', fp);
    }
    return 0;

 fail:
    fprintf(stderr,"%s:%d: Can't write state file\n",
            __FILE__,__LINE__);
    return EIO;
}
