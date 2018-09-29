#include "param.h"
#include "misc.h"
#include "dtnorm.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>

void Param_init(Param *self, const char *name, double value,
                 double low, double high,
                 unsigned type) {
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
/// new itself. Does not allocate "constr" field of
/// constrained parameters. This is set to NULL.
void Param_copy(Param *new, const Param *old) {
    memcpy(new, old, sizeof(Param));
    CHECKMEM(new);
    new->name = strdup(old->name);
    CHECKMEM(new->name);
    if(new->type & CONSTRAINED) {
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
    c = lhs->type - rhs->type;
    if(c)
        return c;
    if(lhs->formula == NULL && rhs->formula==NULL)
        return 0;
    assert(lhs->formula != NULL);
    assert(rhs->formula != NULL);
    return strcmp(lhs->formula, rhs->formula);
}

void Param_sanityCheck(const Param *self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self, file, line);
    REQUIRE(isPow2(self->type & (FREE | FIXED | CONSTRAINED)), file, line);
    unsigned u = self->type & (TWON | TIME | MIXFRAC | ARBITRARY);
    REQUIRE(isPow2(u), file, line);
    REQUIRE(self->name, file, line);
    REQUIRE(legalName(self->name), file, line);
    REQUIRE(isfinite(self->value), file, line);
    REQUIRE(self->low <= self->value, file, line);
    REQUIRE(self->value <= self->high, file, line);
    if(self->type & CONSTRAINED) {
        REQUIRE(self->formula, file, line);
        REQUIRE(self->constr, file, line);
    }else{
        REQUIRE(NULL == self->formula, file, line);
        REQUIRE(NULL == self->constr, file, line);
    }
    if(self->next)
        REQUIRE(self->next > self, file, line);
#endif
}

/// Randomize parameter if FREE and not TIME.
void Param_randomize(Param *self, gsl_rng *rng) {
    assert(self);
    double r;
    if( !(self->type & FREE) )
        return;

    if(self->type & TWON) {
        self->value = dtnorm(self->value, 10000.0, self->low,
                             self->high, rng);
    }else if( self->type & MIXFRAC ) {
        self->value = gsl_ran_beta(rng, 1.0, 5.0);
    }else if( self->type & ARBITRARY ) {
        if(isinf(self->low) && isinf(self->high)) {
            // both bounds are infinite
            self->value += gsl_ran_gaussian_ziggurat(rng, 100.0);
        }else if(isinf(self->low)) {
            // lower bound infinite
            do{
                r = gsl_ran_gaussian_ziggurat(rng, 100.0);
            }while(self->value + r > self->high);
            self->value += r;
        }else if(isinf(self->high)) {
            // upper bound infinite
            do{
                r = gsl_ran_gaussian_ziggurat(rng, 100.0);
            }while(self->value + r < self->low);
            self->value += r;
        }else {
            // finite bounds
            self->value = dtnorm(self->value, 100.0, self->low,
                                 self->high, rng);
        }
    }
}
