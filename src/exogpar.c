#include "exogpar.h"
#include "misc.h"
#include "misc.h"
#include <gsl/gsl_randist.h>

struct ExogPar {
    double *valptr;
    double mean, sd;
};

ExogPar *ExogPar_new(double *valptr, mean, sd) {
    ExogPar *self = malloc(sizeof *self);
    CHECKMEM(self);

    self->valptr = valptr;
    *self->valptr = self->mean = mean;
    self->sd = sd;
    return self;
}

void ExogPar_free(ExogPar *self) {
    free(self);
}

double ExogPar_sample(ExogPar *self, double low, double high, gsl_rng *rng) {
    double x;
    assert(self->mean >= low);
    assert(self->mean <= high);
    if(self->sd == 0.0)
        return *self->valptr;
    do {
        x = gsl->mean + gsl_ran_gaussian(rng, self->sd);
    }while(x <= low || x >= high);
    *self->valptr = x;
    return x;
}


