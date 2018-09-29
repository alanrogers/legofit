#include "param.h"
#include "misc.h"
#include <string.h>
#include <stdlib.h>

void Param_init(Param *self, const char *name, double value,
                 double low, double high,
                 Behavior behavior) {
    assert(self);
    if(low > value || high < value) {
        fprintf(stderr,"%s:%d: can't initialize parameter \"%s\".\n"
                " Value (%g) is not in [%lg, %lg]\n",
                __FILE__,__LINE__, name, value, low, high);
        exit(EXIT_FAILURE);
    }
    self->name = strdup(name);
    CHECKMEM(self->name);
    self->value = value;
    self->low = low;
    self->high = high;
    self->behavior = behavior;
    self->formula = NULL;
    self->constr = NULL;
    self->next = NULL;
}

/// Copy old into new, allocating all the internal pointers but not 
/// new itself. Does not allocate "constr" field of
/// constrained parameters. This is set to NULL.
void Param_copy(Param *new, const Param *old) {
    memcpy(new, old, sizeof(Param));
    CHECKMEM(new);
    new->name = strdup(old->name);
    CHECKMEM(new->name);
    if(new->behavior == Constrained) {
        assert(old->formula != NULL);
        new->formula = strdup(old->formula);
        CHECKMEM(new->formula);
    }else {
        assert(old->formula == NULL);
        new->formula = NULL;
    }
    new->constr = NULL;
    new->next = NULL;
}

// Push a new Param onto the end of a linked list.
Param *Param_push(Param *self, Param *new) {
    if(self == NULL)
        return new;
    self->next = Param_push(self->next, new);
    return self;
}

// frees only memory allocated within Param, not Param itself
void Param_freePtrs(Param *self) {
    free(self->name);
    if(self->formula)
        free(self->formula);
    if(self->constr)
        te_free(self->constr);
}

/// Print name and value of a Param.
void Param_print(Param *self, FILE *fp) {
    fprintf(fp, "   %8s = %lg\n", self->name, self->value);
}

int Param_compare(const Param *lhs, const Param *rhs) {
    int c;
    c = strcmp(lhs->name, rhs->name);
    if(c)
        return c;
    if(lhs->value > rhs->value)
        return 1;
    if(lhs->value < rhs->value)
        return -1;
    if(lhs->low > rhs->low)
        return 1;
    if(lhs->low < rhs->low)
        return -1;
    if(lhs->high > rhs->high)
        return 1;
    if(lhs->high < rhs->high)
        return -1;
    c = lhs->behavior - rhs->behavior;
    if(c)
        return c;
    if(lhs->formula == NULL && rhs->formula==NULL)
        return 0;
    assert(lhs->formula != NULL);
    assert(rhs->formula != NULL);
    return strcmp(lhs->formula, rhs->formula);
}

