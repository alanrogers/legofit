/**
 * @file idsetset.c
 * @author Alan R. Rogers
 * @brief A table of IdSet values.
 *
 * The IdSetSet_add function merges entries that have identical tipId_t values
 * and also identical migration histories. In that case, the probabilities of
 * duplicate IdSet values add. If the IdSet objects are not identical,
 * the IdSetSet_add function as separate copies.
 *
 * To iterate across an IdSetSet named x:
 *
 * IdSet *idset;
 * IdSetSet_rewind(x);
 * while( (idset = IdSetSet_next(x)) != NULL ) {
 *     <operate on idset>
 * }
 *
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "idsetset.h"
#include "idset.h"
#include "error.h"
#include "binary.h"
#include "misc.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/errno.h>

static int resize(IdSetSet *self, int dim);

#define SETS_PER_BUCKET 0.7

static const double sets_per_bucket = SETS_PER_BUCKET;
static const double buckets_per_set = 1.0/SETS_PER_BUCKET;

typedef struct El El;

/// A single key-value pair, with a pointer to the next one.
/// Within each bin of the hash table, key-value pairs are
/// stored in a sorted linked list.
struct El {
    struct El  *next;          ///< next item in list
    IdSet *idset;
};

/// The hash table.
struct IdSetSet {
    int      dim;   // current size of table
    unsigned mask;
    long     nelem; // number of elements stored
    long     maxelem; // max nelem before resize
    El     **tab;   // array of dim pointers to El pointers.
    El *curr;
    int currndx;
};

static El  *El_new(IdSet *idset);
static El  *El_free_shallow(El * e);
static El  *El_free_deep(El * e);
static El  *El_add(El *self, IdSet *idset, int *status);
static void El_sanityCheck(El *el, const char *file, int line);

/// Construct a new element with given key and value
static El *El_new(IdSet *idset) {
    El *new = malloc(sizeof(*new));
    CHECKMEM(new);

    new->next = NULL;
    new->idset = idset;
    return new;
}

/// Add an new IdSet into the linked list. Usage:
///
///     El *list=NULL;
///     int status;
///     list = El_add(list, idset, &status);
///     switch(status) {
///     case 0:
///         printf("Added to probability of existing IdSet.\n");
///         break;
///     case 1:
///         printf("Added new IdSet to list\n");
///         break;
///     default:
///         printf("This should not happen.\n");
///     }
/// @param self current element of list
/// @param[in] idset, a pointer
/// @param[out] status pointer to int, which will be set to 1 if the
/// new IdSet did not previously exist in the list or to 0 if it was
/// already there.
/// @return pointer to self or to a newly-allocated El.
static El *El_add(El *self, IdSet *idset, int *status) {
    El *new;
    if(self == NULL) {
        new = El_new(idset);
        new->next = NULL;
        *status = 1;
        return new;
    }

    int cmp = IdSet_cmp(self->idset, idset);
    if(cmp < 0) {
        new = El_new(idset);
        new->next = self;
        *status = 1;
        return new;
    }else if(cmp > 0) {
        self->next = El_add(self->next, idset, status);
        return self;
    }else{
        // Identical IdSets. If this ever happens, abort.
        fprintf(stderr,"%s:%s:%d: cannot add identical IdSet objects\n",
                __FILE__,__func__,__LINE__);
        IdSet_print(self->idset, stderr);
        putc('\n', stderr);
        IdSet_print(idset, stderr);
        putc('\n', stderr);
        exit(EXIT_FAILURE);
    }
    IdSet_free(idset);
    *status = 0;
    return self;
}

/// Destroy a linked list of El objects, but do not free IdSet pointers.
static El *El_free_shallow(El * e) {
    if(e == NULL)
        return NULL;
    El_free_shallow(e->next);
    free(e);
    return NULL;
}

/// Destroy a linked list of El objects, including IdSet pointers.
static El *El_free_deep(El * e) {
    if(e == NULL)
        return NULL;
    El_free_deep(e->next);
    IdSet_free(e->idset);
    free(e);
    return NULL;
}

static void El_sanityCheck(El *el, const char *file, int line) {
    if(el == NULL)
        return;
    El_sanityCheck(el->next, file, line);
    IdSet_sanityCheck(el->idset, file, line);
}

/// Constructor reserves space for n elements.
IdSetSet *IdSetSet_new(int n) {
    IdSetSet  *new = malloc(sizeof(*new));
    CHECKMEM(new);

    int dim = next_power_of_2(n * buckets_per_set);
    new->dim = dim;
    new->mask = dim - 1;
    new->maxelem = ceil(sets_per_bucket * new->dim);
    new->nelem = 0;
    new->tab = malloc(dim * sizeof(new->tab[0]));
    CHECKMEM(new->tab);

    memset(new->tab, 0, dim * sizeof(new->tab[0]));
    new->curr = NULL;
    new->currndx = -1;
    return new;
}

/// Shallow destructor does not free IdSet objects
void IdSetSet_free_shallow(IdSetSet * self) {
    for(int i=0; i < self->dim; ++i)
        self->tab[i] = El_free_shallow(self->tab[i]);
    free(self->tab);
    free(self);
}

/// Deep destructor frees IdSet objects
void IdSetSet_free_deep(IdSetSet * self) {
    for(int i=0; i < self->dim; ++i)
        self->tab[i] = El_free_deep(self->tab[i]);
    free(self->tab);
    free(self);
}

/// Empty IdSetSet. Does not free the IdSet pointers.
void IdSetSet_empty_shallow(IdSetSet * self) {
    for(int i=0; i < self->dim; ++i)
        self->tab[i] = El_free_shallow(self->tab[i]);
    self->nelem = 0;
    self->curr = NULL;
    self->currndx = -1;
}

/// Empty IdSetSet, including IdSet pointers.
void IdSetSet_empty_deep(IdSetSet * self) {
    for(int i=0; i < self->dim; ++i)
        self->tab[i] = El_free_deep(self->tab[i]);
    self->nelem = 0;
    self->curr = NULL;
    self->currndx = -1;
}

/// Add a value to the table, resizing if necessary.  @return 0 on
/// success; ENOMEM if the function attempts unsuccessfully to resize
/// the hash table.
int IdSetSet_add(IdSetSet * self, IdSet *idset) {

    int status;

    if(1+self->nelem > self->maxelem) {
        status = resize(self, 2 * self->dim);
        if(status)
            return ENOMEM;
    }

    unsigned h = IdSet_hash(idset) & self->mask;
    assert(h < self->dim);

    assert(self);

    self->tab[h] = El_add(self->tab[h], idset, &status);
    self->nelem += status;
    
    return 0;
}

/// Return the number of elements.
int IdSetSet_size(IdSetSet * self) {
    return self->nelem;
}

/// Put idset pointers into array "v". If the array isn't large enough,
/// return BUFFER_OVERFLOW. Otherwise return 0.
int IdSetSet_toArray(IdSetSet *self, unsigned size, IdSet *v[size]) {
    int box, j;
    for(box=j=0; box < self->dim; ++box) {
        El *el;
        for(el = self->tab[box]; el; el = el->next) {
            if(j == size)
                return BUFFER_OVERFLOW;
            v[j++] = el->idset;
        }
    }
    return 0;
}

/// Change the size of the hash map and rehash.
static int resize(IdSetSet *self, int dim) {
    if(!isPow2(dim)) {
        fprintf(stderr,"%s:%d: dim=%d is not a power of 2\n",
                __FILE__,__LINE__, dim);
        exit(EXIT_FAILURE);
    }
    int status;
    unsigned mask = dim-1;
    El **tab = malloc(dim * sizeof(tab[0]));
    if(tab == NULL)
        return ENOMEM;
    memset(tab, 0, dim * sizeof(tab[0]));

    for(int i=0; i < self->dim; ++i) {
        for(El *e=self->tab[i]; e; e=e->next) {
            unsigned long h = IdSet_hash(e->idset) & mask;
            assert(h < dim);
            tab[h] = El_add(tab[h], e->idset, &status);
            switch(status) {
            case 0:
                fprintf(stderr,"%s:%d: duplicate key\n",__FILE__,__LINE__);
                free(tab);
                exit(EXIT_FAILURE);
                break;
            case 1:
                break;
            default:
                fprintf(stderr,"%s:%d: unknown error\n",__FILE__,__LINE__);
                free(tab);
                exit(EXIT_FAILURE);
                break;
            }
        }
        self->tab[i] = El_free_shallow(self->tab[i]);
    }
    free(self->tab);
    self->tab = tab;
    self->dim = dim;
    self->maxelem = ceil(sets_per_bucket * self->dim);
    self->mask = mask;

#ifndef NDEBUG
    IdSetSet_sanityCheck(self,__FILE__,__LINE__);
#endif    
    return 0;
}

/// Move curr to first filled bucket in hash table. Return 0
/// on success or 1 if the hash table is empty;
int IdSetSet_rewind(IdSetSet *self) {
    int i;
    for(i=0; i < self->dim; ++i) {
        if(self->tab[i] != NULL) {
            self->currndx = i;
            self->curr = self->tab[i];
            break;
        }
    }
    if(i == self->dim) {
        self->currndx = -1;
        self->curr = NULL;
        return 1; // failure
    }
    return 0;
}

/// Return the next IdSet pointer, or NULL if there are no more. Then
/// move to next entry in the hash table.
IdSet *IdSetSet_next(IdSetSet *self) {
    if(self->curr == NULL)
        return NULL;
    IdSet *rval = self->curr->idset;
    self->curr = self->curr->next;
    if(self->curr == NULL) {
        // Look for a bucket that isn't empty.
        int i;
        for(i = self->currndx+1; i < self->dim; ++i) {
            if(self->tab[i]) {
                self->currndx = i;
                self->curr = self->tab[i];
                break;
            }
        }

        // All remaining buckets are empty.
        if(i == self->dim) {
            self->currndx = -1;
            self->curr = NULL;
        }
    }
    return rval;
}

/// Make sure hash table is large enough to hold m additional
/// elements.
int IdSetSet_reserve(IdSetSet *self, int m) {
    m += self->nelem;
    if(m > self->maxelem) {
        m = next_power_of_2(m * buckets_per_set);
        int status = resize(self, m);
        if(status)
            return ENOMEM;
    }
    return 0;
}

void IdSetSet_sanityCheck(IdSetSet *self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self != NULL, file, line);
    REQUIRE(self->dim > 0, file, line);
    REQUIRE(isPow2(self->dim), file, line);
    REQUIRE(self->mask == self->dim - 1, file, line);
    REQUIRE(self->nelem >= 0, file, line);
    REQUIRE(self->nelem <= self->maxelem, file, line);
    REQUIRE(self->maxelem == ceil(sets_per_bucket * self->dim), file, line);
    REQUIRE(self->tab != 0, file, line);
    REQUIRE(self->currndx == -1 || (self->currndx >=0 &&
                                    self->currndx < self->dim),
            file, line);
    for(int i=0; i < self->dim; ++i)
        El_sanityCheck(self->tab[i], file, line);
#endif    
}

