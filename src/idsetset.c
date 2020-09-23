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
 * IdSetSet_rewind(x);
 * for( IdSet *idset = IdSetSet_next(x);
 *      idset;
 *      idset = IdSetSet_next(x); {
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

static int resize(IdSetSet *self);

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
    El     **tab;   // array of dim pointers to El pointers.
    El *curr;
    int currndx;
};

static El  *El_new(IdSet *idset);
static void El_free(El * e);
static El  *El_add(El *self, IdSet *idset, int *status);

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
    }
    // Identical IdSets. Add the probability of the new
    // one to that of the old one and free the new one.
    self->idset->p += idset->p;
    IdSet_free(idset);
    *status = 0;
    return self;
}

/// Destroy a linked list of El objects.
static void El_free(El * e) {
    if(e == NULL)
        return;
    El_free(e->next);
    IdSet_free(e->idset);
    free(e);
}

/// Constructor
IdSetSet *IdSetSet_new(int dim) {
    IdSetSet  *new = malloc(sizeof(*new));
    CHECKMEM(new);

    dim = next_power_of_2(dim);
    new->dim = dim;
    new->mask = dim - 1;
    new->nelem = 0;
    new->tab = malloc(dim * sizeof(new->tab[0]));
    if(new->tab == NULL) {
        fprintf(stderr,"%s:%d: bad malloc\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    memset(new->tab, 0, dim * sizeof(new->tab[0]));
    new->curr = NULL;
    new->currndx = -1;
    return new;
}

/// Destructor
void IdSetSet_free(IdSetSet * self) {
    for(int i=0; i < self->dim; ++i)
        El_free(self->tab[i]);
    free(self->tab);
    free(self);
}

/// Empty IdSetSet, freeing all contained IdSet objects,
/// but not freeing the IdSetSet.
void IdSetSet_empty(IdSetSet * self) {
    for(int i=0; i < self->dim; ++i)
        El_free(self->tab[i]);
    self->nelem = 0;
    self->curr = NULL;
    self->currndx = -1;
}

/// Add a value to the table, resizing if necessary.  @return 0 on
/// success; ENOMEM if the function attempts unsuccessfully to resize
/// the hash table.
int IdSetSet_add(IdSetSet * self, IdSet *idset) {

    int status;

    if(self->nelem > 0.7 * self->dim) {
        status = resize(self);
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
/// I may need to make this faster.
int IdSetSet_size(IdSetSet * self) {
    int size = 0;

    for(int i = 0; i < self->dim; ++i) {
        El *el;
        for(el = self->tab[i]; el; el = el->next)
            ++size;
    }
    return size;
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

/// Double the size of the hash map and rehash.
static int resize(IdSetSet *self) {
    int dim = 2*self->dim, status;
    unsigned mask = dim-1;
    El **tab = malloc(dim * sizeof(tab[0]));
    if(tab == NULL)
        return ENOMEM;

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
        El_free(self->tab[i]);
    }
    free(self->tab);
    self->tab = tab;
    self->dim = dim;
    self->mask = mask;
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
