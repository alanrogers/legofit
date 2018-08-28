#include "param.h"
#include "misc.h"
#include <string.h>
#include <stdlib.h>

Param *Param_new(const char *name, double *valptr, double low, double high,
                 ParamStatus status) {
    Param *self = malloc(sizeof(Param));
    CHECKMEM(self);
    self->name = strdup(name);
    self->valptr = valptr;
    self->low = low;
    self->high = high;
    self->status = status;
    return self;
}

void Param_free(Param *self) {
    free(self->name);
    free(self);
}
