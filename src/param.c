#include "param.h"
#include "misc.h"
#include <string.h>
#include <stdlib.h>

Param *Param_new(const char *name, double value, double low, double high,
                 unsigned type) {
    Param *self = malloc(sizeof(Param));
    CHECKMEM(self);
    self->name = strdup(name);
    self->value = value;
    self->low = low;
    self->high = high;
    self->type = type;
    return self;
}

void Param_free(Param *self) {
    free(self->name);
    free(self);
}
