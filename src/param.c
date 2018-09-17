#include "param.h"
#include "misc.h"
#include <string.h>
#include <stdlib.h>

void Param_init(Param *self, const char *name, double value,
                 double low, double high,
                 unsigned type) {
    assert(self);
    if(low > value || high < value) {
        fprintf(stderr,"%s:%d: can't initialize parameter \"%s\"."
                " Value (%g) is not in [%lg, %lg]\n",
                __FILE__,__LINE__, name, low, high);
        exit(EXIT_FAILURE);
    }
    self->name = strdup(name);
    self->value = value;
    self->low = low;
    self->high = high;
    self->type = type;
    self->mean = self->sd = strtod("NaN", NULL);
}

void Param_setGaussianParams(Param *self, double mean, double sd) {
    assert(self);
    self->mean = mean;
    self->sd = sd;
}

// frees only memory allocated within Param, not Param itself
void Param_freePtrs(Param *self) {
    free(self->name);
}

/// Print name and value of a Param if it is of type "onlytype"
void Param_print(Param *self, unsigned onlytype; FILE *fp) {
    if(self && (self->type & onlytype))
        fprintf(fp, "   %8s = %lg\n", self->name, self->value);
}

