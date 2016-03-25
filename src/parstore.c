/**
 * @file parstore.c
 * @brief Manage a vector of parameters
 *
 * This class solves the following problem. One algorithm requires
 * that parameters be distributed throughout a network of
 * nodes. Another requires that they all be collected into a
 * vector. To accomplish both goals, ParStore maintains a vector of
 * parameter values, together with high and low bounds on those
 * values. When a new parameter is added to the ParStore, a pointer to
 * that value is returned, so that it can be stored in the distributed
 * data structure.
 *
 * @copyright Copyright (c) 2016, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "parstore.h"

/// Constructor
ParStore   *ParStore_new(void) {
    ParStore   *self = malloc(sizeof(ParStore));
    checkmem(self, __FILE__, __LINE__);
    self->n = 0;
    return self;
}

/// Destructor
void ParStore_free(ParStore * self) {
    free(self);
}

/// Add parameter to ParStore and return a pointer to that value.
double     *ParStore_addPar(ParStore * self, double value, double lo, double hi) {
    int         i = self->n;

    if(++self->n >= MAXPAR)
        eprintf("%s:%s:%d: buffer overflow in ParStore. "
                "Increase MAXPAR and recompile.\n",
                __FILE__, __func__, __LINE__);

    if(value < lo || value >= hi)
        eprintf("%s:%s:%d: value (%lf) not in range [%lf,%lf)\n",
                __FILE__, __func__, __LINE__, value, lo, hi);

    self->val[i] = value;
    self->lo[i] = lo;
    self->hi[i] = hi;
    return self->val + i;
}
