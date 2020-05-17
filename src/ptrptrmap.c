/**
 * @file ptrptrmap.c
 * @author Alan R. Rogers
 * @brief Map pointers to pointers.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "ptrptrmap.h"
#include "binary.h"
#include "misc.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define KEYSIZE 20

/// PNT_DIM: dimention of hash table must be a power of 2
#define PNT_DIM 32u

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
    const void *key;
    void *value;
};

/// The hash table.
struct PtrPtrMap {
    El         *tab[PNT_DIM];
};

El  *El_new(const void *key, void *value);
void El_free(El * e);
El  *El_insert(El *self, const void *key, void *value, int *status);
static unsigned long ptrHash( const void * key);

static unsigned long ptrHash( const void * key) {
    if(sizeof(void *) == 4)
        return uint32Hash((uint32_t) key);
    else if(sizeof(void *)==8)
        return uint64Hash((uint64_t) key);
    else {
        fprintf(stderr,"%s:%s:%d: can't hash pointer\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }
    return 0;
}

/// Construct a new element with given key and value
El *El_new(const void *key, void *value) {
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
El *El_insert(El *self, const void *key, void *value, int *status) {
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

/// PtrPtrMap constructor
PtrPtrMap    *PtrPtrMap_new(void) {
    PtrPtrMap    *new = malloc(sizeof(*new));
    CHECKMEM(new);
    memset(new->tab, 0, sizeof(new->tab));
    return new;
}

/// PtrPtrMap destructor
void PtrPtrMap_free(PtrPtrMap * self) {
    int i;
    for(i=0; i < PNT_DIM; ++i)
        El_free(self->tab[i]);
    free(self);
}

/// Get the value associated with key. If no such object exists,
/// return NULL.
void *PtrPtrMap_get(PtrPtrMap * self, const void *key) {

    // Same as hash % PNT_DIM but faster. Requires
    // that PNT_DIM be a power of 2.
    unsigned long h = ptrHash(key) & (PNT_DIM-1ul);
    assert(h < PNT_DIM);

    assert(self);

    El *el;
    for(el = self->tab[h]; el; el = el->next) {
        if(key == el->key)
            return el->value;
        else if(key < el->key)
            return NULL;
    }
    return NULL;
}

/// Insert a pointer into the table.
/// @return 0 on success; 1 if pointers is already in table.
int PtrPtrMap_insert(PtrPtrMap * self, const void *key, void *value) {

    // Same as hash % PNT_DIM but faster. Requires
    // that PNT_DIM be a power of 2.
    unsigned long h = ptrHash(key) & (PNT_DIM-1ul);
    assert(h < PNT_DIM);
    int status;

    assert(self);

    self->tab[h] = El_insert(self->tab[h], key, value, &status);
    return status;
}

/// Return the number of elements in the PtrPtrMap.
unsigned long PtrPtrMap_size(PtrPtrMap * self) {
    unsigned    i;
    unsigned long size = 0;

    for(i = 0; i < PNT_DIM; ++i) {
        El *el;
        for(el = self->tab[i]; el; el = el->next)
            ++size;
    }
    return size;
}

/// Print a PtrPtrMap
void PtrPtrMap_print(PtrPtrMap *self) {
    unsigned i;
    for(i=0; i < PNT_DIM; ++i) {
        printf("%2u:", i);
        El *el;
        for(el = self->tab[i]; el; el = el->next)
            printf(" [%p => %p]", el->key, el->value);
        putchar('\n');
    }
}

