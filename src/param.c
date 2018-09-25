#include "param.h"
#include "misc.h"
#include <string.h>
#include <stdlib.h>

void Param_init(Param *self, const char *name, double value,
                 double low, double high,
                 ParamType type) {
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
    self->type = type;
    self->formula = NULL;
    self->constr = NULL;
    self->next = NULL;
}

/// Copy old into new, allocating all the internal pointers but not 
/// new itself.
void Param_copy(Param *new, const Param *old, const te_variables *variables,
                int *error) {
    memcpy(new, old, sizeof(Param));
    CHECKMEM(new);
    new->name = strcpy(old->name);
    CHECKMEM(new->name);
    if(new->type == Constrained) {
        new->formula = strcpy(old->formula);
        CHECKMEM(new->formula);
        new->constr = te_compile(new->formula, variables, error);
        CHECKMEM(new->constr);
    }else{
        new->formula = NULL;
        new->constr = NULL;
    }
    new->next = NULL;
}

// Push a new Param onto the end of a linked list.
Param *Param_push(Param *self, Param *new) {
    if(self == NULL)
        return new;
    fprintf(stderr,"%s:%s:%d\n", __FILE__,__func__,__LINE__);
    self->next = Param_push(self->next, new);
    fprintf(stderr,"%s:%s:%d\n", __FILE__,__func__,__LINE__);
    return self;
}

// frees only memory allocated within Param, not Param itself
void Param_freePtrs(Param *self) {
    fprintf(stderr,"%s:%s:%d\n", __FILE__,__func__,__LINE__);
    free(self->name);
    if(self->formula) {
        fprintf(stderr,"%s:%s:%d\n", __FILE__,__func__,__LINE__);
        free(self->formula);
    }
    if(self->constr) {
        fprintf(stderr,"%s:%s:%d\n", __FILE__,__func__,__LINE__);
        te_free(self->constr);
    }
}

/// Print name and value of a Param if it is of type "onlytype"
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
    c = lhs->type - rhs->type;
    if(c)
        return c;
    if(lhs->formula == NULL && rhs->formula==NULL)
        return 0;
    assert(lhs->formula != NULL);
    assert(rhs->formula != NULL);
    return strcmp(lhs->formula, rhs->formula);
}

