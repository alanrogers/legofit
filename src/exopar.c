#include "exopar.h"
#include "misc.h"
#include "misc.h"
#include <gsl/gsl_randist.h>

struct ExoPar {
    double *valptr;
    double mean, sd;
};

ExoPar *ExoPar_new(double *valptr, mean, sd) {
    ExoPar *self = malloc(sizeof *self);
    CHECKMEM(self);

    self->valptr = valptr;
    *self->valptr = self->mean = mean;
    self->sd = sd;
    return self;
}

void ExoPar_free(ExoPar *self) {
    free(self);
}

double ExoPar_sample(ExoPar *self, double low, double high, gsl_rng *rng) {
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


