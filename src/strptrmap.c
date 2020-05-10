/**
 * @file strptrmap.c
 * @author Alan R. Rogers
 * @brief Map strings to pointers.
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "strptrmap.h"
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
    struct El  *next;          ///< ptr to next item in list
    char        key[KEYSIZE];  
    void       *ptr;
};

/// The hash table.
struct StrPtrMap {
    El         *tab[PNT_DIM];
};

El         *El_new(const char *key, void *ptr);
void        El_free(El * e);
El         *El_insert(El *self, const char *key, void *ptr, int *status);

/// Construct a new element with given key and pointer
El         *El_new(const char *key, void *ptr) {
    El         *new = malloc(sizeof(*new));
    CHECKMEM(new);

    new->next = NULL;
    snprintf(new->key, sizeof(new->key), "%s", key);
    assert(0 == strncmp(new->key, key, sizeof(new->key)));
    new->ptr = ptr;
    return new;
}

/// Insert a new key/ptr pair into the linked list. Usage:
///
///     El *list=NULL;
///     int status;
///     list = El_insert(list, key, ptr, &status);
///     if(status != 0)
///         printf("Error: key is already in list\n");
///
/// @param self current element of list
/// @param[in] key, a character string
/// @param[in] ptr pointer 
/// @param[out] status pointer to int, which will be set to 0 on
/// success or 1 if ptr is already in list.
/// @return pointer to self or to a newly-allocated El.
El *El_insert(El *self, const char *key, void *ptr, int *status) {
    El *new;
    if(self == NULL) {
        new = El_new(key, ptr);
        new->next = NULL;
        *status = 0; // success
        return new;
    }
    int cmp = strncmp(key, self->key, sizeof(new->key));
    if(cmp < 0) {
        new = El_new(key, ptr);
        new->next = self;
        *status = 0; // success
        return new;
    }else if(cmp > 0) {
        self->next = El_insert(self->next, key, ptr, status);
        return self;
    }
    assert(cmp==0);
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

/// StrPtrMap constructor
StrPtrMap    *StrPtrMap_new(void) {
    StrPtrMap    *new = malloc(sizeof(*new));
    CHECKMEM(new);
    memset(new->tab, 0, sizeof(new->tab));
    return new;
}

/// StrPtrMap destructor
void StrPtrMap_free(StrPtrMap * self) {
    int i;
    for(i=0; i < PNT_DIM; ++i)
        El_free(self->tab[i]);
    free(self);
}

/// Get the pointer associated with key. If no such object exists,
/// return NULL.
void *StrPtrMap_get(StrPtrMap * self, const char *key) {

    // Same as hash % PNT_DIM but faster. Requires
    // that PNT_DIM be a power of 2.
    unsigned long h = strhash(key) & (PNT_DIM-1ul);
    assert(h < PNT_DIM);

    assert(self);

    El *el;
    for(el = self->tab[h]; el; el = el->next) {
        int cmp = strncmp(key, el->key, sizeof(el->key));
        if(cmp == 0)
            return el->ptr;
        else if(cmp < 0)
            return NULL;
    }
    return NULL;
}

/// Insert a pointer into the table.
/// @return 0 on success; 1 if pointers is already in table.
int StrPtrMap_insert(StrPtrMap * self, const char *key,
                      void *ptr) {

    // Same as hash % PNT_DIM but faster. Requires
    // that PNT_DIM be a power of 2.
    unsigned long h = strhash(key) & (PNT_DIM-1ul);
    assert(h < PNT_DIM);
    int status;

    assert(self);

    self->tab[h] = El_insert(self->tab[h], key, ptr, &status);
    return status;
}

/// Return the number of elements in the StrPtrMap.
unsigned long StrPtrMap_size(StrPtrMap * self) {
    unsigned    i;
    unsigned long size = 0;

    for(i = 0; i < PNT_DIM; ++i) {
        El *el;
        for(el = self->tab[i]; el; el = el->next)
            ++size;
    }
    return size;
}

/// Print a StrPtrMap
void StrPtrMap_print(StrPtrMap *self) {
    unsigned i;
    for(i=0; i < PNT_DIM; ++i) {
        printf("%2u:", i);
        El *el;
        for(el = self->tab[i]; el; el = el->next)
            printf(" [%s, %p]", el->key, el->ptr);
        putchar('\n');
    }
}

/// Fill array v with pointers. Abort if v is of wrong size.
void StrPtrMap_ptrArray(StrPtrMap *self, long unsigned n, void *v[n]) {
    unsigned j=0;
    for(unsigned i=0; i < PNT_DIM; ++i) {
        El *el;
        for(el = self->tab[i]; el; el = el->next) {
            if(j == n) {
                fprintf(stderr,"%s:%s:%d: output array is of size %lu;"
                        " should %lu\n",
                        __FILE__,__func__,__LINE__, n, StrPtrMap_size(self));
                exit(EXIT_FAILURE);
            }
            v[j++] = el->ptr;
        }
    }
    if(j != n)  {
        fprintf(stderr,"%s:%s:%d: output array is of size %lu;"
                " should %lu\n", 
                __FILE__,__func__,__LINE__, n, StrPtrMap_size(self));
        exit(EXIT_FAILURE);
    }
}


