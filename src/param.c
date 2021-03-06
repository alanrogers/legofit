#include "param.h"
#include "misc.h"
#include "dtnorm.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>

int Param_isFree(const Param *self) {
    return self->type & FREE;
}

double Param_getValue(const Param *self) {
    return self->value;
}

/// If value is within the range allowed for this parameter,
/// then set the Param to the provided value and return 0.
/// Otherwise leave the Param unchanged and return EDOM.
int Param_setValue(Param *self, double value) {
    if(value < self->low || value > self->high)
        return EDOM;

    self->value = value;
    return 0;
}

Param *Param_new(const char *name, double value,
               double low, double high,
               unsigned type, const char *formula) {
    if(low > value || high < value) {
        fprintf(stderr,"%s:%d: can't initialize parameter \"%s\".\n"
                " Value (%g) is not in [%lg, %lg]\n",
                __FILE__,__LINE__, name, value, low, high);
        exit(EXIT_FAILURE);
    }
    Param *self = malloc(sizeof(Param));
    CHECKMEM(self);

    assert(name);
    self->name = strdup(name);
    CHECKMEM(self->name);
    self->value = value;
    self->low = low;
    self->high = high;
    self->type = type;
    if(formula)
        self->formula = strdup(formula);
    else
        self->formula = NULL;
    self->constr = NULL;
    return self;
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
    }else{
        REQUIRE(NULL == self->formula, file, line);
    }
#endif
}

/// If parameter is FREE return a random value within the
/// range of legal values. Otherwise, return self->value.
double Param_getTrialValue(const Param *self, gsl_rng *rng) {
    assert(self);
    double trial, stdev;
    if( !(self->type & FREE) )
        return self->value;

    // trial value
    if(self->type & TWON) {
        stdev = self->value == 0 ? 1.0 : 2.0 * self->value;
        trial = dtnorm(self->value, stdev, self->low, self->high, rng);
    }else if( self->type & TIME ) {
        if(isfinite(self->low)) {
            stdev = self->value == 0 ? 1.0 : 2.0 * self->value;
            trial = dtnorm(self->value, stdev, self->low, self->high, rng);
        }else {
            assert( !isfinite(self->low) );
            fprintf(stderr, "%s:%s:%d: non-finite lower bound of TIME"
                    " parameter\n", __FILE__,__func__,__LINE__);
            exit(EXIT_FAILURE);
        }
    }else if( self->type & MIXFRAC ) {
        trial = gsl_ran_beta(rng, 1.0, 5.0);
    }else if( self->type & ARBITRARY ) {
        if(isinf(self->low) && isinf(self->high)) {
            // both bounds are infinite
            trial = self->value + gsl_ran_gaussian_ziggurat(rng, 100.0);
        }else if(isinf(self->low)) {
            // lower bound infinite
            do{
                trial = gsl_ran_gaussian_ziggurat(rng, 100.0);
            }while(self->value + trial > self->high);
            trial += self->value;
        }else if(isinf(self->high)) {
            // upper bound infinite
            do{
                trial = gsl_ran_gaussian_ziggurat(rng, 100.0);
            }while(self->value + trial < self->low);
            trial += self->value;
        }else {
            // finite bounds
            stdev = 0.5*(self->high - self->low);
            trial = dtnorm(self->value, stdev, self->low,
                                 self->high, rng);
        }
    }else{
        fprintf(stderr,"%s:%d: invalid Param type %o\n",
                __FILE__,__LINE__,self->type);
        exit(EXIT_FAILURE);
    }

    return trial;
}

/// Move the contents of "from" to "to".
void Param_move(Param *to, Param *from) {
    assert(to);
    assert(from);
    Param_freePtrs(to);
    memcpy(to, from, sizeof(Param));
    memset(from, 0, sizeof(Param));
}

/// Copy from into to, but don't copy from->constr, which
/// must be set by a separate call to Param_compileConstraint.
void Param_copy(Param *to, const Param *from) {
    assert(to);
    assert(from);
    Param_freePtrs(to);
    memcpy(to, from, sizeof(Param));
    to->name = strdup(from->name);
    CHECKMEM(to->name);
    if(from->formula) {
        to->formula = strdup(from->formula);
        CHECKMEM(to->formula);
    }
    to->constr = NULL;
}

void Param_compileConstraint(Param *self, StrPtrMap *te_pars) {

    int status;

    // Ignore parameters that aren't constrained
    if( (self->type & CONSTRAINED) == 0)
        return;
    assert(self->formula);
    
    // compile constraint
    self->constr = te_compile(self->formula, te_pars, &status);
    if(self->constr == NULL) {
        fprintf(stderr, "%s:%d: parse error\n", __FILE__, __LINE__);
        fprintf(stderr, "  %s\n", self->formula);
        fprintf(stderr, "  %*s^\nError near here\n", status - 1, "");
        exit(EXIT_FAILURE);
    }
}

/// Set value of constrained parameter.  Abort with an error message
/// if result is NaN.
void Param_constrain(Param *par) {
#ifndef NDEBUG    
    if(par->constr == NULL) {
        fprintf(stderr,"%s:%s:%d: NULL argument %p\n",
                __FILE__,__func__,__LINE__, par->constr);
        fprintf(stderr,"formula: %s = %s\n",
                par->name, par->formula);
        exit(EXIT_FAILURE);
    }
#endif    
    par->value = te_eval(par->constr);
    if(isnan(par->value)) {
        fprintf(stderr,"%s:%d: constraint returned NaN\n",
                __FILE__,__LINE__);
        fprintf(stderr,"formula: %s = %s\n",
                par->name, par->formula);
        exit(EXIT_FAILURE);
    }
}

// frees only memory allocated within Param, not Param itself
void Param_freePtrs(Param *self) {
    assert(self);
    if(self->name)
        free(self->name);
    if(self->formula)
        free(self->formula);
    if(self->constr)
        te_free(self->constr);
    self->name = self->formula = NULL;
    self->constr = NULL;
}
