#ifndef ARR_PARSTORE_H
# define ARR_PARSTORE_H

#  include "typedefs.h"
#  include "misc.h"
#  define MAXPAR 10

struct ParStore {
    int n;                // number of parameters
    double val[MAXPAR];   // parameter values
    double lo[MAXPAR];    // lower bounds
    double hi[MAXPAR];    // upper bounds
};

ParStore *ParStore_new(void);
void ParStore_free(ParStore *self);
double *addPar(ParStore *self, double value, double lo, double hi);

ParStore *ParStore_new(void) {
    ParStore *self = malloc(sizeof(ParStore));
    checkmem(self, __FILE__, __LINE__);
    self->n = 0;
    return self;
}

void ParStore_free(ParStore *self) {
    free(self);
}

/// Add parameter to ParStore and return a pointer to that value.
double *addPar(ParStore *self, double value, double lo, double hi) {
    int i = self->n;

    if(++self->n >= MAXPAR) {
        eprintf("%s:%s:%d: buffer overflow in ParStore. "
                "Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__);
    }

    self->val[i] = value;
    self->lo[i] = lo;
    self->hi[i] = hi;
    return self->val + i;
}

#endif
