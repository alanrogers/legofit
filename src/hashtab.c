/**
 * @file hashtab.c
 * @author Alan R. Rogers
 * @brief Hash table mapping keys (character strings) to values (void
 * pointers).
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "hashtab.h"
#include "misc.h"
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define KEYSIZE 20

/// HASHDIM: dimention of hash table must be a power of 2
#define HASHDIM 32u

// Make sure HASHDIM is a power of 2
#if (HASHDIM==0u || (HASHDIM & (HASHDIM-1u)))
# error HASHDIM must be a power of 2
#endif

/// A single key-value pair, with a pointer to the next one.
/// Within each bin of the hash table, key-value pairs are
/// stored in a sorted linked list.
struct El {
    struct El  *next;
    char        key[KEYSIZE];
    void       *value;
};

/// The hash table.
struct HashTab {
    El         *tab[HASHDIM];
};

/// Treats the hash table as a continuous sequence.
struct HashTabSeq {
    struct HashTab *ht;
    int ndx;
    El *el;
};

El         *El_new(const char *key);
void        El_free(El * e);
El         *El_find(El * self, El ** found, const char *key);
void        El_map(El *self, void(*f)(void *value, void *data), void *data);
void        HashTabSeq_nextBucket(HashTabSeq *self);
void        freeValue(void *value, void *unused);

/// For use by HashTab_map. Will free every non-NULL value in table.
void freeValue(void *p, void *unused) {
    if(p)
        free(p);
}

/// Construct a new element with given key and NULL value.
El         *El_new(const char *key) {
    El         *new = malloc(sizeof(*new));
    CHECKMEM(new);

    new->next = NULL;
    snprintf(new->key, sizeof(new->key), "%s", key);
    assert(0 == strncmp(new->key, key, sizeof(new->key)));
    new->value = NULL;
    return new;
}

/// Destroy a linked list of El objects.
void El_free(El * e) {
    if(e == NULL)
        return;
    El_free(e->next);
    free(e);
}

/// Return the value pointer of an El object.
void * El_get(El * self) {
    assert(self);
    return self->value;
}

/// Set the value pointer of an El object
void El_set(El * self, void *value) {
    assert(self);
    self->value = value;
}

/// Return El corresponding to key. If no such El exists,
/// allocate a new one, insert it into hashtab, and return that.
/// Note that El_find does not return NULL if key is missing.
El *El_find(El * self, El ** found, const char *key) {
    int comparison;
    if(self == NULL
       || 0 > (comparison = strncmp(key, self->key, sizeof(self->key)))) {
        *found = El_new(key);
        (*found)->next = self;
        return *found;
    } else if(comparison > 0) {
        self->next = El_find(self->next, found, key);
        return self;
    }
    assert(0 == strncmp(key, self->key, sizeof(self->key)));
    *found = self;
    return self;
}

/// Print an El but not the entire linked list.
void El_printShallow(El * self) {
    if(self == NULL)
        return;
    printf(" [%s, %p]", self->key, self->value);
}

/// Print a linked list of El objects.
void El_print(El * self) {
    if(self == NULL)
        return;
    El_printShallow(self);
    El_print(self->next);
}

/// Apply function f(value, data) to all values in list.
void El_map(El *self, void(*f)(void *value, void *data), void *data) {
    if(self==NULL)
        return;
    (*f)(self->value, data);
    El_map(self->next, f, data);
}

/// HashTab constructor
HashTab    *HashTab_new(void) {
    HashTab    *new = malloc(sizeof(*new));
    CHECKMEM(new);
    memset(new->tab, 0, sizeof(new->tab));
    return new;
}

/// HashTab destructor
void HashTab_free(HashTab * self) {
    int i;
    for(i=0; i < HASHDIM; ++i)
        El_free(self->tab[i]);
    free(self);
}

/// Get the El object associated with key. If no
/// such object exists in HashTab, insert a new one with
/// NULL value pointer and return that.
El *HashTab_get(HashTab * self, const char *key) {

    // Same as hash % HASHDIM but faster. Requires
    // that HASHDIM be a power of 2.
    unsigned h = strhash(key) & (HASHDIM-1u);
    assert(h < HASHDIM);

    assert(self);

    El *el = NULL;
    self->tab[h] = El_find(self->tab[h], &el, key);
    assert(el);

    return el;
}

/// Return the number of elements in the HashTab.
unsigned long HashTab_size(HashTab * self) {
    unsigned    i;
    unsigned long size = 0;

    for(i = 0; i < HASHDIM; ++i) {
        El *el;
        for(el = self->tab[i]; el; el = el->next)
            ++size;
    }
    return size;
}

/// Print a HashTab
void HashTab_print(HashTab *self) {
    unsigned i;
    for(i=0; i < HASHDIM; ++i) {
        printf("%2u:", i);
        El_print(self->tab[i]);
        putchar('\n');
    }
}

/// Apply function f(value, data) to all values in table.
void HashTab_map(HashTab *self, void(*f)(void *value, void *data), void *data) {
    unsigned i;

    for(i=0; i<HASHDIM; ++i)
        El_map(self->tab[i], f, data);
}

/// HashTabSeq constructor.
HashTabSeq *HashTabSeq_new(HashTab *ht) {
    HashTabSeq *self = malloc(sizeof(*self));
    CHECKMEM(self);

    self->ht = ht;
    self->ndx = -1;
    self->el = NULL;
    return self;
}

/// Find the next bucket, ignoring empty buckets. If all remaining buckets
/// are empty, set self->ndx equal to HASHDIM.
void HashTabSeq_nextBucket(HashTabSeq *self) {
    ++self->ndx;
    while(self->ndx < HASHDIM && self->ht->tab[self->ndx]==NULL)
        ++self->ndx;
}

/// Return the next element in table.
El *HashTabSeq_next(HashTabSeq *self) {
    if(self->ndx < 0) {
        HashTabSeq_nextBucket(self);
        if(self->ndx == HASHDIM) {
            self->el = NULL;
            return NULL;
        }
        self->el = self->ht->tab[self->ndx];
        return self->el;
    }
    if(self->el == NULL) {
        if(self->ndx != HASHDIM) {
            fprintf(stderr,"%s:%s:%d: this shouldn't happen.",
                    __FILE__, __func__, __LINE__);
            exit(EXIT_FAILURE);
        }
        return NULL;
    }
    self->el = self->el->next;
    if(self->el == NULL) {
        HashTabSeq_nextBucket(self);
        if(self->ndx == HASHDIM)
            return NULL;
        self->el = self->ht->tab[self->ndx];
    }
    return self->el;
}

/// HashTabSeq destructor
void HashTabSeq_free(HashTabSeq *self) {
    free(self);
}
