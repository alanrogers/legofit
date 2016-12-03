#include "exopar.h"
#include "misc.h"
#include <string.h>
#include <assert.h>
#include <gsl/gsl_randist.h>

struct ExoPar {
    double *ptr; // not locally owned
    double mean, sd, low, high;
};

struct ExoParTab {
    int n;
    ExoPar *v;
};

struct ExoParList {
    double *ptr;  // not locally owned
    ExoPar par;
};

void ExoPar_init(ExoPar *self, double *ptr, double mean, double sd,
                 double low, double high);
int ExoParList_size(ExoParList *self);
int compare_ExoPar_ExoPar(const void *void_x, const void void_y);
int compare_dblPtr_ExoPar(const void *void_x, const void void_y);

void ExoPar_init(ExoPar *self, double *ptr, double mean, double sd,
                 double low, double high) {
    self->ptr = ptr;
    *self->ptr = self->mean = mean;
    self->sd = sd;
    self->low = low;
    self->high = high;
}

double ExoPar_sample(ExoPar *self, double low, double high, gsl_rng *rng) {
    double x;
    assert(self->mean >= low);
    assert(self->mean <= high);
    if(self->sd == 0.0)
        return *self->ptr;
    do {
        x = gsl->mean + gsl_ran_gaussian(rng, self->sd);
    }while(x <= low || x >= high);
    *self->ptr = x;
    return x;
}

/// Add a new link to list.
ExoParList *ExoParList_add(ExoParList *old, double *ptr, double m, double sd,
                           double sd, double low, double high) {
    ExoParList *new = malloc(sizeof ExoParList);
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

int compare_ExoPar_ExoPar(const void *void_x, const void void_y) {
    const ExoPar * x = (const ExoPar *) void_x;
    const ExoPar * y = (const ExoPar *) void_y;

    if(x->ptr > y->ptr)
        return 1;
    else if(x->ptr < y->ptr)
        return -1;
    assert(x->ptr==y->ptr);
    return 0;
}

int compare_dblPtr_ExoPar(const void *void_x, const void void_y) {
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
ExoParTab *ExoParTab_new(ExoParList *list) {
    ExoParTab *self = malloc(sizeof ExoParTab);
    CHECKMEM(self);

    self->n = ExoParList_size(list);
    self->v = malloc(self->n * sizeof(ExoPar));
    CHECKMEM(self->v);
    for(int i=0; i < self->n; ++i) {
        assert(list != NULL);
        memcpy(self->v + i, &list->par, sizeof ExoPar);
        list = list->next;
    }
    qsort(self->v, (size_t) self->n, sizeof(ExoPar),
          compare_ExoPar_ExoPar);
    return self;
}

/// Return const pointer to ExoPar object corresponding to
/// ptr, or NULL if no such object is found.
ExoPar const * const ExoParTab_find(ExoParTab *self, double *ptr) {
    return bsearch(&ptr, self->v, (size_t) self->n,
                   sizeof ExoPar, compare_dblPtr_ExoPar);
}
