/**
 * @file idset_tbl.c
 * @author Alan R. Rogers
 * @brief A table of IdSet values.
 *
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "idset_tbl.h"
#include "error.h"
#include "binary.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/errno.h>

static int resize(IdSetTbl *self);

typedef struct El El;

/// A single key-value pair, with a pointer to the next one.
/// Within each bin of the hash table, key-value pairs are
/// stored in a sorted linked list.
struct El {
    struct El  *next;          ///< next item in list
    IdSet *idset;
};

/// The hash table.
struct IdSetTbl {
    int      dim;   // current size of table
    unsigned mask;
    long     nelem; // number of elements stored
    El     **tab;   // array of dim pointers to El pointers.
};

static El  *El_new(IdSet *idset)
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
///         printf("Added new IdSet to list\n");
///         break;
///     case 1:
///         printf("Added to probability of existing IdSet.\n");
///         break;
///     default:
///         printf("This should not happen.\n");
///     }
/// @param self current element of list
/// @param[in] idset, a pointer
/// @param[out] status pointer to int, which will be set to 0 if the
/// new IdSet did not previously exist in the list or to 1 if it was
/// already there.
/// @return pointer to self or to a newly-allocated El.
static El *El_add(El *self, IdSet *idset, int *status) {
    El *new;
    if(self == NULL) {
        new = El_new(idset);
        new->next = NULL;
        *status = 0;
        return new;
    }

    int cmp = IdSet_cmp(self->idset, idset);
    if(cmp < 0) {
        new = El_new(idset);
        new->next = self;
        *status = 0;
        return new;
    }else if(cmp > 0) {
        self->next = El_add(self->next, idset, status);
        return self;
    }
    // Identical IdSets. Add the probability of the new
    // one to that of the old one and free the new one.
    self->idset->p += idset->p;
    IdSet_free(idset);
    *status = 1;
    return self;
}

/// Destroy a linked list of El objects.
static void El_free(El * e) {
    if(e == NULL)
        return;
    El_free(e->next);
    free(e);
}

/// Constructor
IdSetTbl *IdSetTbl_new(int dim) {
    IdSetTbl  *new = malloc(sizeof(*new));
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
    return new;
}

/// Destructor
void IdSetTbl_free(IdSetTbl * self) {
    for(int i=0; i < self->dim; ++i)
        El_free(self->tab[i]);
    free(self->tab);
    free(self);
}

/// Add a value to the table, resizing if necessary.
/// @return 0 on success; 1 if pointers is already in table; ENOMEM
/// if the function attempts unsuccessfully to resize the hash table.
int IdSetTbl_add(IdSetTbl * self, IdSet *idset) {

    int status;

    if(self->nelem > 0.7 * self->dim) {
        status = resize(self);
        if(status)
            return status;
    }

    unsigneed h = IdSet_hash(idset) & self->mask;
    assert(h < self->dim);

    assert(self);

    self->tab[h] = El_add(self->tab[h], idset, &status);
    if(status == 1)
        self->nelem += 1;
    
    return status;
}

/// Return the number of elements.
unsigned long FUNC(IdSetTbl, size)(IdSetTbl * self) {
    unsigned    i;
    unsigned long size = 0;

    for(i = 0; i < self->dim; ++i) {
        El *el;
        for(el = self->tab[i]; el; el = el->next)
            ++size;
    }
    return size;
}

/// Put idsets into array "idsets". If the array isn't large enough,
/// return BUFFER_OVERFLOW. Otherwise return 0.
int FUNC(IdSetTbl, keys)(IdSetTbl *self, unsigned size, IdSet *keys[size]) {
    unsigned box, j;
    for(box=j=0; box < self->dim; ++box) {
        El *el;
        for(el = self->tab[box]; el; el = el->next) {
            if(j == size)
                return BUFFER_OVERFLOW;
            keys[j++] = el->key;
        }
    }
    return 0;
}

/// Double the size of the hash map and rehash.
static int resize(IdSetTbl *self) {
    int dim = 2*self->dim, status;
    unsigned mask = dim-1;
    El **tab = malloc(dim * sizeof(tab[0]));
    if(tab == NULL)
        return ENOMEM;

    for(int i=0; i < self->dim; ++i) {
        for(El *e=self->tab[i]; e; e=e->next) {
            unsigned long h = HASH( (KEYCAST) e->idset) & mask;
            assert(h < dim);
            tab[h] = El_add(tab[h], e->idset, e->value, &status);
            if(status)
                goto error;
        }
        El_free(self->tab[i]);
    }
    free(self->tab);
    self->tab = tab;
    self->dim = dim;
    self->mask = mask;
    return 0;
    
 error:
    free(tab);
    return status;
}
