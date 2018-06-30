#include "state.h"
#include "gptree.h"
#include "error.h"
#include "misc.h"
#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

typedef enum StateFile_t StateFile_t;
enum StateFile_t {UNSET, OLD, NEW};

struct NameList {
    char *name;
    struct NameList *next;
};

struct State {
    StateFile_t filetype;
    int npts, npar; // numbers of points and parameters
    int haveNames;  // boolean: nonzero means names have been assigned
    char **name;     // name[j] is name of j'th parameter
    double *cost;   // cost[i] is cost function at i'th point
    double **s;     // s[i][j]=value of j'th param at i'th point
};

NameList *NameList_append(NameList *self, const char *name) {
    if(self == NULL) {
        self = malloc(sizeof(NameList));
        CHECKMEM(self);
        self->name = strdup(name);
        self->next = NULL;
    }else {
        if(0 == strcmp(name, self->name)) {
            fprintf(stderr, "%s:%d: state file \"%s\" listed multiple times.\n",
                    __FILE__,__LINE__,name);
            exit(EXIT_FAILURE);
        }
        self->next = NameList_append(self->next, name);
    }
    return self;
}

void NameList_free(NameList *self) {
    if(self == NULL)
        return;
    NameList_free(self->next);
    free(self);
}

int NameList_size(NameList *self) {
    if(self==0)
        return 0;
    return 1 + NameList_size(self->next);
}

void NameList_print(NameList *self, FILE *fp) {
    if(self==NULL)
        return;
    fprintf(fp, " %s", self->name);
    NameList_print(self->next, fp);
}

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
    self->name = malloc(npar * sizeof(self->name[0]));
    CHECKMEM(self->name);
    memset(self->name, 0, npar * sizeof(self->name[0]));
    self->haveNames = 0;
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
    for(i=0; i < self->npar; ++i) {
        if(self->name[i])
            free(self->name[i]);
    }
    free(self->name);
    free(self->s);
    free(self);
}

// Construct a new State object by reading a file
State *State_read(FILE *fp) {
    int i, j, npts, npar, status;
    State *self = NULL;

    {
        // Read first line and figure out whether we're in
        // a new-format state file or an old-format file.
        char buff[200], fmt[200], dummy[200];
        if(fgets(buff, sizeof(buff), fp) == NULL) {
            fprintf(stderr,"%s:%d: empty state file\n",__FILE__,__LINE__);
            goto fail;
        }
        if(NULL == strchr(buff, '\n') && !feof(fp)) {
            fprintf(stderr,"%s:%d: Buffer overflow. size=%zu\n",
                    __FILE__,__LINE__, sizeof(buff));
            goto fail;
        }
        status = sscanf(buff, "%d %d %s %s", &npts, &npar, fmt, dummy);
        switch(status) {
        case 2:
            self->filetype = OLD;
            break;
        case 3:
            self->filetype = NEW;
            if(0 != strcmp("new_format", fmt)) {
                fprintf(stderr,"%s:%d: state file format error\n",
                        __FILE__,__LINE__);
                fprintf(stderr,"  Old format: 1st line has 2 ints\n");
                fprintf(stderr,"  New format: 2 ints, then \"new_format\"\n");
                fprintf(stderr,"  Got \"%s\" instead of \"new_format\"\n",
                        fmt);
                goto fail;
            }
            break;
        default:
            fprintf(stderr,"%s:%d: state file format error\n",
                    __FILE__,__LINE__);
            fprintf(stderr,"  Old format: 1st line has 2 fields\n");
            fprintf(stderr,"  New format: 3 fields\n");
            fprintf(stderr,"  Input: %s", buff);
            goto fail;
        }
    }
    
    self = State_new(npts, npar);
    CHECKMEM(self);

    if(self->filetype == NEW) {
        // Read header containing parameter names
        char buff[100];
        status = fscanf(fp, "%s", buff);
        if(status != 1) {
            fprintf(stderr,"%s:%d: status=%d\n", __FILE__,__LINE__,status);
            goto fail;
        }
        if(0 != strcmp("cost", buff)) {
            fprintf(stderr, "%s:%d: got \"%s\" instead of \"cost\""
                    " reading new-format state file\n",
                    __FILE__,__LINE__, buff);
            goto fail;
        }
        for(j=0; j < npar; ++j) {
            status = fscanf(fp, "%s", buff);
            if(status != 1) {
                fprintf(stderr,"%s:%d: status=%d\n", __FILE__,__LINE__,status);
                goto fail;
            }
            if(!isalpha(buff[0])) {
                fprintf(stderr,"%s:%d: parameter names must start with"
                        " an alphabetic character\n", __FILE__,__LINE__);
                fprintf(stderr,"  Got \"%s\"\n", buff);
                goto fail;
            }
            self->name[j] = strdup(buff);
            if(self->name[j] == NULL) {
                fprintf(stderr,"%s:%d: bad strdup\n", __FILE__,__LINE__);
                goto fail;
            }
        }
    }

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

/// Set name of parameter with index "ndx". Return 0 on success,
/// abort on failure.
int State_setName(State *self, int ndx, const char *name) {
    assert(ndx < self->npar);
    assert(ndx >= 0);

    self->name[ndx] = strdup(name);
    if(self->name[ndx] == NULL) {
        fprintf(stderr,"%s:%d: bad strdup\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    return 0;
}

// Print State object to a file
int State_print(State *self, FILE *fp) {
    int i, j, status, imin=0;

    // find index of minimum cost
    for(i=1; i < self->npts; ++i)
        if(self->cost[i] < self->cost[imin])
            imin = i;

    status = fprintf(fp, "%d %d new_format\n", self->npts, self->npar);
    if(status==0)
        goto fail;

    // write header containing column names
    assert(self->haveNames);
    status = fprintf(fp, "%s", "cost");
    if(status==0)
        goto fail;
    for(j=0; j < self->npar; ++j) {
        status = fprintf(fp, " %s", self->name[j]);
        if(status == 0)
            goto fail;
    }
    status = putc('\n', fp);
    if(status == EOF)
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

    // print remaining lines
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
State *State_readList(NameList *list, int npts, GPTree *gptree) {
    int nstates = NameList_size(list);
    if(nstates==0)
        return NULL;
    int npar = GPTree_nFree(gptree);

    State *state[nstates];
    NameList *node;
    int i, j, k;

    for(i=0, node=list; i<nstates; ++i, node=node->next) {

        // Create a State object from each file name in list.
        FILE *fp = fopen(node->name, "r");
        state[i] = State_read(fp);
        CHECKMEM(state[i]);
        fclose(fp);

        // Make sure dimensions are compatible.
        if(npar != State_nparameters(state[i])) {
            fprintf(stderr,"%s:%s:%d:"
                    " input state file \"%s\" has"
                    " incompatible dimensions.\n",
                    __FILE__,__func__,__LINE__,
                    node->name);
            exit(EXIT_FAILURE);
        }

        if(state[i]->filetype != NEW) {
            // Can't check parameter names in old-format
            // state files.
            continue;
        }

        // Make sure parameter names are consistent
        for(j=0; j < npar; ++j) {
            if(0 != strcmp(GPTree_getNameFree(gptree, j),
                           state[i]->name[j])) {
                fprintf(stderr,"%s:%s:%d:"
                        " input state file \"%s\" has"
                        " incompatible parameter names.\n",
                        __FILE__,__func__,__LINE__,
                        node->name);
                fprintf(stderr,"    \"%s\" != \"%s\"\n",
                        GPTree_getNameFree(gptree, j),
                        state[i]->name[j]);
                exit(EXIT_FAILURE);
            }
        }
    }

    // Make sure we're not asking for more points than exist
    int total=0;
    for(i=0; i < nstates; ++i)
        total += State_npoints(state[i]);
    if(total < npts)
        npts = total;

    State *self = State_new(npts, npar);
    CHECKMEM(self);

    // Figure out how many points to take from each State object.
    // The goal is to take nearly equal numbers.
    div_t qr = div(npts, nstates);
    int npoints[nstates];
    int avail[nstates];
    int got=0;

    // 1st pass: roughly equal number of points per input file
    for(i=0; i < nstates; ++i) {
        avail[i] = State_npoints(state[i]);
        npoints[i] = qr.quot;
        if(i < qr.rem)
            ++npoints[i];

        if(npoints[i] > avail[i])
            npoints[i] = avail[i];
        got += npoints[i];
    }

    // additional passes: fill in missing points
    while(got < npts) {
        for(i=0; i < nstates; ++i) {
            if(npoints[i] < avail[i]) {
                ++npoints[i];
                ++got;
            }
            if(got == npts)
                break;
        }
    }

    // Set vectors in self; free pointers in vector "state".
    k = 0; // index into self->s
    for(i=0; i < nstates; ++i) {
        for(j=0; j < npoints[i]; ++j) {
            State_setCost(self, k, state[i]->cost[j]);
            State_setVector(self, k, npar, state[i]->s[j]);
            ++k;
        }
        State_free(state[i]);
    }
    return self;
}
