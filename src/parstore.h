#ifndef ARR_PARSTORE_H
#  define ARR_PARSTORE_H

#  include "typedefs.h"
#  include "misc.h"
#  include <assert.h>
#  define MAXPAR 10

struct ParStore {
    int         n;              // number of parameters
    double      val[MAXPAR];    // parameter values
    double      lo[MAXPAR];     // lower bounds
    double      hi[MAXPAR];     // upper bounds
};

struct Bounds {
	double lo_twoN, hi_twoN, lo_t, hi_t;
};

ParStore   *ParStore_new(void);
void        ParStore_free(ParStore * self);
double     *ParStore_addPar(ParStore * self, double value, double lo, double hi);

static inline int ParStore_nPar(ParStore * self);
static inline double ParStore_get(ParStore * self, int i);
static inline void ParStore_set(ParStore * self, int i, double value);
static inline double ParStore_loBnd(ParStore * self, int i);
static inline double ParStore_hiBnd(ParStore * self, int i);
static inline double *ParStore_getPtr(ParStore * self);

/// Return the number of parameters
static inline int ParStore_nPar(ParStore * self) {
    return self->n;
}

/// Get value of i'th parameter
static inline double ParStore_get(ParStore * self, int i) {
    assert(i < self->n);

    return self->val[i];
}

/// Set value of i'th parameter
static inline void ParStore_set(ParStore * self, int i, double value) {
    assert(i < self->n);

    self->val[i] = value;
}

/// Return low bound of i'th parameter
static inline double ParStore_loBnd(ParStore * self, int i) {
    assert(i < self->n);

    return self->lo[i];
}

/// Return high bound of i'th parameter
static inline double ParStore_hiBnd(ParStore * self, int i) {
    assert(i < self->n);

    return self->hi[i];
}

/// Return pointer to array of values
static inline double * ParStore_getPtr(ParStore * self) {
    return &self->val[0];
}

#endif
