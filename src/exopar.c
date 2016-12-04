#include "exopar.h"
#include "misc.h"
#include <string.h>
#include <assert.h>
#include <gsl/gsl_randist.h>

struct ExoPar {
    double *ptr; // not locally owned
    double mean, sd, low, high;
};

struct ExoParList {
    double *ptr;  // not locally owned
    ExoPar par;
    ExoParList *next;
};

struct ExoParTab {
    ExoParList *list;
    int frozen;
    int n;
    ExoPar *v;
};

static void ExoPar_init(ExoPar *self, double *ptr, double mean, double sd,
                        double low, double high);
static double ExoPar_sample(const ExoPar *self, double low, double high,
                            gsl_rng *rng);
int ExoParList_size(ExoParList *self);
ExoParList *ExoParList_add(ExoParList *old, double *ptr, double m, double sd,
                           double low, double high);
static int compare_ExoPar_ExoPar(const void *void_x, const void *void_y);
static int compare_dblPtr_ExoPar(const void *void_x, const void *void_y);

static void ExoPar_init(ExoPar *self, double *ptr, double mean, double sd,
                 double low, double high) {
    self->ptr = ptr;
    *self->ptr = self->mean = mean;
    self->sd = sd;
    self->low = low;
    self->high = high;
}

static double ExoPar_sample(const ExoPar *self, double low, double high,
                            gsl_rng *rng) {
    double x;
    assert(self->mean >= low);
    assert(self->mean <= high);
    if(self->sd == 0.0)
        return *self->ptr;
    do {
        x = self->mean + gsl_ran_gaussian(rng, self->sd);
    }while(x <= low || x >= high);
    *self->ptr = x;
    return x;
}

/// Add a new link to list.
// ptr points to the memory occupied by the parameter; m is the mean,
// sd the standard deviation, low the lower bound, and high the upper
// bound.
ExoParList *ExoParList_add(ExoParList *old, double *ptr, double m, 
                           double sd, double low, double high) {
    ExoParList *new = malloc(sizeof(ExoParList));
    CHECKMEM(new);

    ExoPar_init(&new->par, ptr, m, sd, low, high);
    new->next = old;
    return new;
}

void ExoParList_free(ExoParList *self) {
    if(self==NULL)
        return;
    ExoParList_free(self->next);
    free(self);
}

/// Count links in list
int ExoParList_size(ExoParList *self) {
    int n=0;
    while(self != NULL) {
        ++n;
        self = self->next;
    }
    return n;
}

static int compare_ExoPar_ExoPar(const void *void_x, const void *void_y) {
    const ExoPar * x = (const ExoPar *) void_x;
    const ExoPar * y = (const ExoPar *) void_y;

    if(x->ptr > y->ptr)
        return 1;
    else if(x->ptr < y->ptr)
        return -1;
    assert(x->ptr==y->ptr);
    return 0;
}

static int compare_dblPtr_ExoPar(const void *void_x, const void *void_y) {
    double * const * x = (double * const *) void_x;
    const ExoPar * y = (const ExoPar *) void_y;

    if(*x > y->ptr)
        return 1;
    else if(*x < y->ptr)
        return -1;
    assert(*x == y->ptr);
    return 0;
}

/// Create a new ExoParTab object by coping ExoPar values
/// from linked list and then sorting them.
ExoParTab *ExoParTab_new(void) {
    ExoParTab *self = malloc(sizeof(ExoParTab));
    CHECKMEM(self);

    self->list = NULL;
    self->frozen = 0;
    self->n = 0;
    self->v = NULL;
    return self;
}

void ExoParTab_freeze(ExoParTab *self) {

    if(self->frozen)
        DIE("Can't freeze an ExoParTab twice");
    self->frozen = 1;
    self->n = ExoParList_size(self->list);
    self->v = malloc(self->n * sizeof(ExoPar));
    CHECKMEM(self->v);
    for(int i=0; i < self->n; ++i) {
        assert(self->list != NULL);
        memcpy(self->v + i, &self->list->par, sizeof(ExoPar));
        self->list = self->list->next;
    }
    qsort(self->v, (size_t) self->n, sizeof(ExoPar),
          compare_ExoPar_ExoPar);
    ExoParList_free(self->list);
    self->list = NULL;
}

/// Return a new value sampled from the distribution associated with
/// ptr.
double ExoParTab_sample(ExoParTab *self, double *ptr,
                                double low, double high,
                                gsl_rng *rng) {
    assert(self->frozen);
    const ExoPar *exopar = bsearch(&ptr, self->v,
                                   (size_t) self->n,
                                   sizeof(ExoPar),
                                   compare_dblPtr_ExoPar);
    return ExoPar_sample(exopar, low, high, rng);
}

void ExoParTab_free(ExoParTab *self) {
    if(self->frozen) {
        free(self->v);
        assert(self->list == NULL);
    }else{
        ExoParList_free(self->list);
        assert(self->v == NULL);
    }
    free(self);
}

// Add an exogeneous parameter to the table.
// ptr points to the memory occupied by the parameter; m is the mean,
// sd the standard deviation, low the lower bound, and high the upper
// bound.
void ExoParTab_add(ExoParTab *self, double *ptr, double m, double sd,
                   double low, double high) {
    self->list = ExoParList_add(self->list, ptr, m, sd, low, high);
}
