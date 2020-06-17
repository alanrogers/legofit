/**
 * @file uint64uint64map.c
 * @author Alan R. Rogers
 * @brief Map uint64_t to uint64_t
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "u64u64map.h"
#include "binary.h"
#include "misc.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/// PNT_DIM: dimention of hash table must be a power of 2
#define PNT_DIM 1024u

// Make sure PNT_DIM is a power of 2
#if (PNT_DIM==0u || (PNT_DIM & (PNT_DIM-1u)))
# error PNT_DIM must be a power of 2
#endif

typedef struct El El;

/// A single key-value pair, with a pointer to the next one.
/// Within each bin of the hash table, key-value pairs are
/// stored in a sorted linked list.
struct El {
    struct El  *next;          ///< next item in list
    uint64_t key;
    uint64_t value;
};

/// The hash table.
struct U64U64Map {
    El         *tab[PNT_DIM];
};

El  *El_new(const uint64_t key, uint64_t value);
void El_free(El * e);
El  *El_insert(El *self, uint64_t key, uint64_t value, int *status);

/// Construct a new element with given key and value
El *El_new(const uint64_t key, uint64_t value) {
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
El *El_insert(El *self, uint64_t key, uint64_t value, int *status) {
    El *new;
    if(self == NULL) {
        new = El_new(key, value);
        new->next = NULL;
        *status = 0; // success
        return new;
    }

    if(key < self->key) {
        new = El_new(key, value);
        new->next = self;
        *status = 0; // success
        return new;
    }else if(key > self->key) {
        self->next = El_insert(self->next, key, value, status);
        return self;
    }
    assert(key == self->key);
    *status = 1; // failure
    return self;
}

/// Destroy a linked list of El objects.
void El_free(El * e) {
    if(e == NULL)
        return;
    El_free(e->next);
    free(e);
}

/// U64U64Map constructor
U64U64Map    *U64U64Map_new(void) {
    U64U64Map    *new = malloc(sizeof(*new));
    CHECKMEM(new);
    memset(new->tab, 0, sizeof(new->tab));
    return new;
}

/// U64U64Map destructor
void U64U64Map_free(U64U64Map * self) {
    int i;
    for(i=0; i < PNT_DIM; ++i)
        El_free(self->tab[i]);
    free(self);
}

/// Get the value associated with key. Return 0 on success, 1 if
/// there is no such key
int U64U64Map_get(U64U64Map * self, uint64_t key, uint64_t *value) {

    // Same as hash % PNT_DIM but faster. Requires
    // that PNT_DIM be a power of 2.
    unsigned long h = uint64Hash(key) & (PNT_DIM-1ul);
    assert(h < PNT_DIM);

    assert(self);

    El *el;
    for(el = self->tab[h]; el; el = el->next) {
        if(key == el->key) {
            *value = el->value;
            return 0;
        }else if(key < el->key)
            return 1;
    }
    return 1;
}

/// Insert a value into the table.
/// @return 0 on success; 1 if pointers is already in table.
int U64U64Map_insert(U64U64Map * self, uint64_t key, uint64_t value) {

    // Same as hash % PNT_DIM but faster. Requires
    // that PNT_DIM be a power of 2.
    unsigned long h = uint64Hash(key) & (PNT_DIM-1ul);
    assert(h < PNT_DIM);
    int status;

    assert(self);

    self->tab[h] = El_insert(self->tab[h], key, value, &status);
    return status;
}

/// Return the number of elements in the U64U64Map.
unsigned long U64U64Map_size(U64U64Map * self) {
    unsigned    i;
    unsigned long size = 0;

    for(i = 0; i < PNT_DIM; ++i) {
        El *el;
        for(el = self->tab[i]; el; el = el->next)
            ++size;
    }
    return size;
}

/// Print a U64U64Map
void U64U64Map_print(U64U64Map *self, FILE *fp) {
    unsigned i;
    for(i=0; i < PNT_DIM; ++i) {
        printf("%2u:", i);
        El *el;
        for(el = self->tab[i]; el; el = el->next)
            printf(" [%llu => %llu]", el->key, el->value);
        putchar('\n');
    }
}

