/**
 * @brief Exogeneous parameters sampled from truncated normal distribution.
 * @file exopar.c
 * @author Alan R. Rogers
 */
#include "exopar.h"
#include "misc.h"
#include "dtnorm.h"
#include <string.h>
#include <assert.h>
#include <gsl/gsl_randist.h>

typedef struct ExoParList ExoParList;
typedef struct ExoParItem ExoParItem;

/// ExoParItem represents a parameter whose value will be
/// repeatedly changed by sampling from a truncated Gaussian
/// distribution.
struct ExoParItem {
    double     *ptr;       ///< Pointer to variable. Not locally owned
    double      mean;      ///< Mean of underlying Gaussian dist.
    double      sd;        ///< Std deviation of underlying Gaussian.
};

/// A linked list of ExoParItem objects
struct ExoParList {
    ExoParItem  par;
    ExoParList *next;
};

/// A container of ExoParItem objects.
/// Initially, the container isn't frozen, and the items
/// are maintained in a linked list, which makes it easy to add
/// new items to the container. Later, the container can be
/// frozen. After this, it is no longer possible to add things, and
/// the ExoParItem objects are copied into a sorted array. The frozen
/// array behaves like a dictionary, in which the keys are pointers
/// and the values are ExoParItem objects.
struct ExoPar {
    ExoParList *list;
    int         frozen;
    int         n;
    ExoParItem *v;
};

static void ExoParItem_init(ExoParItem * self, double *ptr, double mean,
                            double sd);
static void ExoParItem_sample(const ExoParItem * self, double low,
                              double high, gsl_rng * rng);
static int  ExoParList_size(ExoParList * self);
static ExoParList *ExoParList_add(ExoParList * old, double *ptr, double m,
                                  double sd);
static void ExoParList_free(ExoParList * self);
static int  compare_ExoParItem_ExoParItem(const void *void_x,
                                          const void *void_y);
static int  compare_dblPtr_ExoParItem(const void *void_x, const void *void_y);

/// Initialize an ExoParItem object.
/// @param [inout] ptr pointer to memory that will be controlled by
/// this ExoParItem.
/// @param[in] mean mean of underlying Gaussian distribution
/// @param[in] sd standard deviation of underlying Gaussian
/// @sideeffect Sets *ptr equal to mean
static void ExoParItem_init(ExoParItem * self, double *ptr, double mean,
                            double sd) {
    assert(ptr);
    self->ptr = ptr;
    *self->ptr = self->mean = mean;
    self->sd = sd;
}

/// Set *ptr by sampling from a truncated Gaussian distribution
/// @param[inout] self ExoParItem object to be modified
/// @param[in] low low end of truncation interval
/// @param[in] high high end of truncation interval
/// @param[inout] rng GSL random number generator
static void ExoParItem_sample(const ExoParItem * self, double low,
                              double high, gsl_rng * rng) {
    double      x;
    assert(low < high);

    // Doubly-truncated normal random variate
    x = dtnorm(self->mean, self->sd, low, high, rng);

    *self->ptr = x;
}

/// Add a new link to list.
/// @param[inout] old pointer to head of linked list
/// @param[in] ptr points to the memory occupied by the parameter to be
/// controlled
/// @param[in] m the mean,
/// @param[in] sd the standard deviation.
/// @return pointer to new head of list
static ExoParList *ExoParList_add(ExoParList * old, double *ptr, double m,
                                  double sd) {
    assert(ptr);
    ExoParList *new = malloc(sizeof(ExoParList));
    CHECKMEM(new);

    ExoParItem_init(&new->par, ptr, m, sd);
    new->next = old;
    return new;
}

/// Free linked list
static void ExoParList_free(ExoParList * self) {
    if(self == NULL)
        return;
    ExoParList_free(self->next);
    free(self);
}

/// Count links in list
static int ExoParList_size(ExoParList * self) {
    int         n = 0;
    while(self != NULL) {
        ++n;
        self = self->next;
    }
    return n;
}

/// Compare two objects of type ExoParItem.
/// @return -1, 0, or 1, to indicate that *x < *y, *x == *y, or *x >
/// *y.
static int compare_ExoParItem_ExoParItem(const void *void_x,
                                         const void *void_y) {
    const ExoParItem *x = (const ExoParItem *) void_x;
    const ExoParItem *y = (const ExoParItem *) void_y;

    if(x->ptr > y->ptr)
        return 1;
    else if(x->ptr < y->ptr)
        return -1;
    assert(x->ptr == y->ptr);
    return 0;
}

/// Compare a double pointer to the "ptr" variable within an
/// ExoParItem.
/// @return -1, 0, or 1 to indicate that the double ptr is "<", "==",
/// or ">" the pointer within the ExoParItem.
static int compare_dblPtr_ExoParItem(const void *void_x, const void *void_y) {
    double     *const *x = (double *const *) void_x;
    const ExoParItem *y = (const ExoParItem *) void_y;

    if(*x > y->ptr)
        return 1;
    else if(*x < y->ptr)
        return -1;
    assert(*x == y->ptr);
    return 0;
}

/// Create an empty ExoPar container.
/// @return pointer to newly-allocated ExoPar.
ExoPar     *ExoPar_new(void) {
    ExoPar     *self = malloc(sizeof(ExoPar));
    CHECKMEM(self);

    self->list = NULL;
    self->frozen = 0;
    self->n = 0;
    self->v = NULL;
    return self;
}

/// Freeze an ExoPar. This converts the internal linked list into
/// a an array, which is sorted according to the internal pointer
/// values. After freezing, nothing can be added to the ExoPar.
void ExoPar_freeze(ExoPar * self) {

    if(self->frozen)
        DIE("Can't freeze an ExoPar twice");
    self->frozen = 1;
    self->n = ExoParList_size(self->list);
    self->v = malloc(self->n * sizeof(ExoParItem));
    CHECKMEM(self->v);
    for(int i = 0; i < self->n; ++i) {
        assert(self->list != NULL);
        memcpy(self->v + i, &self->list->par, sizeof(ExoParItem));
        self->list = self->list->next;
    }
    qsort(self->v, (size_t) self->n, sizeof(ExoParItem),
          compare_ExoParItem_ExoParItem);
    ExoParList_free(self->list);
    self->list = NULL;
}

/// If ptr is in ExoPar, then reset the value it points to with a
/// random value generated by sampling from the truncated normal
/// distribution.
/// @param self pointer to ExoPar
/// @param ptr memory address of the variable we will to modify
/// @param low low end of truncation interval
/// @param high high end of truncation interval
/// @param rng GSL random number generator
/// @return 0 on success or 1 if ptr is not in ExoPar.
int ExoPar_sample(ExoPar * self, double *ptr,
                  double low, double high, gsl_rng * rng) {
    assert(self->frozen);
    const ExoParItem *exopar = bsearch(&ptr, self->v,
                                       (size_t) self->n,
                                       sizeof(ExoParItem),
                                       compare_dblPtr_ExoParItem);
    if(exopar) {
        ExoParItem_sample(exopar, low, high, rng);
        return 0;
    }
    return 1;
}

/// ExoPar destructor
void ExoPar_free(ExoPar * self) {
    if(self->frozen) {
        free(self->v);
        assert(self->list == NULL);
    } else {
        ExoParList_free(self->list);
        assert(self->v == NULL);
    }
    free(self);
}

/// Add an exogeneous parameter to the table.
/// ptr points to the memory occupied by the parameter; m is the mean,
/// sd the standard deviation.
void ExoPar_add(ExoPar * self, double *ptr, double m, double sd) {
    assert(ptr);
    self->list = ExoParList_add(self->list, ptr, m, sd);
}

/// Duplicate an ExoPar object.
ExoPar     *ExoPar_dup(const ExoPar * old) {
    if(!old->frozen)
        DIE("Can't dup an unfrozen ExoPar");

    ExoPar     *new = malloc(sizeof(ExoPar));
    CHECKMEM(new);

    new->n = old->n;
    new->frozen = old->frozen;
    new->list = NULL;
    new->v = memdup(old->v, old->n * sizeof(old->v[0]));
    CHECKMEM(new->v);
    return new;
}

/// Shift all pointers within an ExoPar
/// @param[in] offset magnitude of shift
/// @param[in] sign direction of shift
void ExoPar_shiftPtrs(ExoPar * self, size_t offset, int sign) {
    if(!self->frozen)
        DIE("Can't shift pointers in an unfrozen ExoPar");
    int         i;
    for(i = 0; i < self->n; ++i)
        SHIFT_PTR(self->v[i].ptr, offset, sign);
}
