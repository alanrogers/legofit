/**
 * @file strstrmap.c
 * @author Alan R. Rogers
 * @brief Map char* to double
 * @copyright Copyright (c) 2023, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "strstrmap.h"
#include "error.h"
#include "binary.h"
#include "misc.h"
#include "error.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <sys/errno.h>

static int resize(StrStrMap *self);

typedef struct El El;

/// A single key-value pair, with a pointer to the next one.
/// Within each bin of the hash table, key-value pairs are
/// stored in a sorted linked list.
struct El {
    struct El  *next;          ///< next item in list
    char * key;                // locally owned
    char * value;              // locally owned
};

/// The hash table.
struct StrStrMap {
    int      dim;   // current size of table
    unsigned mask;
    long     nelem; // number of elements stored
    El     **tab;   // array of dim pointers to El pointers.
};

static El  *El_new(const char * key, const char * value);
static void El_free(El * e);
static El  *El_insert(El *self, const char * key, const char * value,
                      int *status);

/// Construct a new element with given key and value
static El *El_new(const char * key, const char * value) {
    El         *new = malloc(sizeof(*new));
    CHECKMEM(new);

    new->next = NULL;
    new->key = strdup(key);
    new->value = strdup(value);
    return new;
}

/// Insert a new key/value pair into the linked list. Usage:
///
///     El *list=NULL;
///     int status;
///     list = El_insert(list, key, value, &status);
///     if(status != 0)
///         printf("Error: key is already in list\n");
///
/// @param self current element of list
/// @param[in] key, a pointer
/// @param[in] value, another pointer 
/// @param[out] status pointer to int, which will be set to 0 on
/// success or 1 if value is already in list.
/// @return pointer to self or to a newly-allocated El.
static El *El_insert(El *self, const char * key, const char * value,
                     int *status) {
    El *new;
    if(self == NULL) {
        new = El_new(key, value);
        new->next = NULL;
        *status = 0; // success
        return new;
    }

    int cmp = strcmp(self->key, key);
    if(cmp < 0) {
        new = El_new(key, value);
        new->next = self;
        *status = 0; // success
        return new;
    }else if(cmp > 0) {
        self->next = El_insert(self->next, key, value, status);
        return self;
    }
    assert( 0 == strcmp(key, self->key) );
    *status = 1; // failure
    return self;
}

/// Destroy a linked list of El objects.
static void El_free(El * e) {
    if(e == NULL)
        return;
    free(e->key);
    free(e->value);
    El_free(e->next);
    free(e);
}

/// Constructor
StrStrMap    * StrStrMap_new(int dim) {
    StrStrMap    *new = malloc(sizeof(*new));
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
void StrStrMap_free(StrStrMap * self) {
    for(int i=0; i < self->dim; ++i)
        El_free(self->tab[i]);
    free(self->tab);
    free(self);
}

/// Get the value associated with key. On success, *status is set
/// equal to 0 and the function returns a const pointer to the value
/// associated with the key. If the key is not found, *status is set
/// equal to 1, and the function returns NULL;
const char * StrStrMap_get(StrStrMap * self, const char * key, int *status) {

    unsigned long h = strhash(key) & self->mask;
    assert(h < self->dim);

    assert(self);

    El *el;
    for(el = self->tab[h]; el; el = el->next) {
        int cmp = strcmp(el->key, key);
        if(cmp == 0) {
            *status = 0;
            return el->value;
        }else if(cmp < 0) {
            *status = 1;
            return NULL;
        }
    }
    *status = 1;
    return NULL;
}

/// Return 1 if key is present in map; 0 otherwise.
int StrStrMap_hasKey(StrStrMap *self, const char *key) {
    unsigned long h = strhash(key) & self->mask;
    assert(h < self->dim);

    assert(self);

    El *el;
    for(el = self->tab[h]; el; el = el->next) {
        int cmp = strcmp(el->key, key);
        if(cmp == 0) {
            return 1;
        }else if(cmp < 0) {
            return 0;
        }
    }
    return 0;
}

/// Insert a value into the table, resizing if necessary.
/// @return 0 on success; 1 if pointers is already in table; ENOMEM
/// if the function attempts unsuccessfully to resize the hash table.
int StrStrMap_insert(StrStrMap * self, const char * key, const char * value) {

    int status;
    
    assert(self);

    if(self->nelem > 0.7 * self->dim) {
        status = resize(self);
        if(status)
            return status;
    }

    unsigned long h = strhash(key) & self->mask;
    assert(h < self->dim);


    self->tab[h] = El_insert(self->tab[h], key, value, &status);

    if(status == 0)
        self->nelem += 1;
    
    return status;
}

/// Return the number of elements.
unsigned long StrStrMap_size(StrStrMap * self) {
    unsigned long size = 0;

    for(unsigned i = 0; i < self->dim; ++i) {
        for(El *el = self->tab[i]; el; el = el->next)
            ++size;
    }
    return size;
}

/// Put keys into array "keys". If the array isn't large enough,
/// return BUFFER_OVERFLOW. Otherwise return 0. Allocates memory
/// for each key, so each key should be freed by the calling
/// function.
int StrStrMap_keys(StrStrMap *self, unsigned size, char * keys[size]) {
    unsigned box, j;
    for(box=j=0; box < self->dim; ++box) {
        for(El *el = self->tab[box]; el; el = el->next) {
            if(j == size) {
                // Array "keys" is too small. Free the strings we've
                // already allocated there and return BUFFER_OVERFLOW.
                for(int i=0; i < size; ++i) {
                    free(keys[i]);
                    keys[i] = NULL;
                }
                return BUFFER_OVERFLOW;
            }
            keys[j++] = strdup(el->key);
        }
    }
    return 0;
}

/// Double the size of the hash map and rehash.
static int resize(StrStrMap *self) {
    int dim = 2*self->dim, status;
    unsigned mask = dim-1;
    El **tab = malloc(dim * sizeof(tab[0]));
    if(tab == NULL)
        return ENOMEM;
    memset(tab, 0, dim * sizeof(tab[0]));

    for(int i=0; i < self->dim; ++i) {
        for(El *e=self->tab[i]; e; e=e->next) {
            unsigned long h = strhash(e->key) & mask;
            assert(h < dim);
            tab[h] = El_insert(tab[h], e->key, e->value, &status);
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
