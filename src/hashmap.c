/**
 * @file hashmap.c
 * @author Alan R. Rogers
 * @brief Generic hashmap code.
 *
 * Generic code for a hash map. To create a concrete version, include
 * this into a .c file that defines MAPTYPE, KEYTYPE, VALTYPE, CMP,
 * HASH, and HASH_DIM.
 *
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#define INNER(CLASS, NAME) CLASS ## _ ## NAME
#define FUNC(CLASS, NAME) INNER(CLASS, NAME)

#include "error.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define XSTR(x) STR(x)
#define STR(x) #x

// Make sure HASH_DIM is a power of 2
#if (HASH_DIM==0u || (HASH_DIM & (HASH_DIM-1u)))
#pragma message "Map type:" XSTR(MAPTYPE)
#pragma message "HASH_DIM:" XSTR(HASH_DIM)
# error $MAPTYPE : HASH_DIM must be a power of 2
#endif

typedef struct El El;

/// A single key-value pair, with a pointer to the next one.
/// Within each bin of the hash table, key-value pairs are
/// stored in a sorted linked list.
struct El {
    struct El  *next;          ///< next item in list
    KEYTYPE key;
    VALTYPE value;
};

/// The hash table.
struct MAPTYPE {
    El         *tab[HASH_DIM];
};

static El  *El_new(KEYTYPE key, VALTYPE value);
static void El_free(El * e);
static El  *El_insert(El *self, KEYTYPE key, VALTYPE value, int *status);

/// Construct a new element with given key and value
static El *El_new(KEYTYPE key, VALTYPE value) {
    El         *new = malloc(sizeof(*new));
    CHECKMEM(new);

    new->next = NULL;
    new->key = key;
    new->value = value;
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
static El *El_insert(El *self, KEYTYPE key, VALTYPE value, int *status) {
    El *new;
    if(self == NULL) {
        new = El_new(key, value);
        new->next = NULL;
        *status = 0; // success
        return new;
    }

    int cmp = CMP(self->key, key);
    if(cmp < 0) {
        new = El_new(key, value);
        new->next = self;
        *status = 0; // success
        return new;
    }else if(cmp > 0) {
        self->next = El_insert(self->next, key, value, status);
        return self;
    }
    assert(key == self->key);
    *status = 1; // failure
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
MAPTYPE    * FUNC(MAPTYPE, new)(void) {
    MAPTYPE    *new = malloc(sizeof(*new));
    CHECKMEM(new);
    memset(new->tab, 0, sizeof(new->tab));
    return new;
}

/// Destructor
void FUNC(MAPTYPE, free)(MAPTYPE * self) {
    int i;
    for(i=0; i < HASH_DIM; ++i)
        El_free(self->tab[i]);
    free(self);
}

/// Get the value associated with key. On success, *status is set equal to 0
/// and the function returns the value associated with the key. If the key
/// is not found, *status is set equal to 1, and the function returns 0,
/// which C will then try to convert into VALTYPE. 
VALTYPE FUNC(MAPTYPE, get)(MAPTYPE * self, const KEYTYPE key, int *status) {

    // Same as hash % HASH_DIM but faster. Requires
    // that HASH_DIM be a power of 2.
    unsigned long h = HASH(key) & (HASH_DIM-1ul);
    assert(h < HASH_DIM);

    assert(self);

    El *el;
    for(el = self->tab[h]; el; el = el->next) {
        int cmp = CMP(el->key, key);
        if(cmp == 0) {
            *status = 0;
            return el->value;
        }else if(cmp < 0) {
            *status = 1;
            return 0;
        }
    }
    *status = 1;
    return 0;
}

/// Insert a value into the table.
/// @return 0 on success; 1 if pointers is already in table.
int FUNC(MAPTYPE, insert)(MAPTYPE * self, KEYTYPE key, VALTYPE value) {

    // Same as hash % HASH_DIM but faster. Requires
    // that HASH_DIM be a power of 2.
    unsigned long h = HASH(key) & (HASH_DIM-1ul);
    assert(h < HASH_DIM);
    int status;

    assert(self);

    self->tab[h] = El_insert(self->tab[h], key, value, &status);
    return status;
}

/// Return the number of elements.
unsigned long FUNC(MAPTYPE, size)(MAPTYPE * self) {
    unsigned    i;
    unsigned long size = 0;

    for(i = 0; i < HASH_DIM; ++i) {
        El *el;
        for(el = self->tab[i]; el; el = el->next)
            ++size;
    }
    return size;
}

/// Put keys into array "keys". If the array isn't large enough,
/// return BUFFER_OVERFLOW. Otherwise return 0.
int FUNC(MAPTYPE, keys)(MAPTYPE *self, unsigned size, KEYTYPE keys[size]) {
    unsigned box, j;
    for(box=j=0; box < HASH_DIM; ++box) {
        El *el;
        for(el = self->tab[box]; el; el = el->next) {
            if(j == size)
                return BUFFER_OVERFLOW;
            keys[j++] = el->key;
        }
    }
    return 0;
}

