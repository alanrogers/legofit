/**
 * @file popnodetab.c
 * @author Alan R. Rogers
 * @brief Hash table associating names of PopNode objects with pointers to
 * PopNode objects.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "popnodetab.h"
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

/// A single key-value pair, with a pointer to the next one.
/// Within each bin of the hash table, key-value pairs are
/// stored in a sorted linked list.
struct El {
    struct El  *next;          ///< ptr to next item in list
    char        key[KEYSIZE];  ///< name of PopNode
    PopNode    *node;          ///< ptr to PopNode; not locally owned
};

/// The hash table.
struct PopNodeTab {
    El         *tab[PNT_DIM];
};

El         *El_new(const char *key, PopNode *node);
void        El_free(El * e);
El         *El_insert(El *self, const char *key, PopNode *node, int *status);

/// Construct a new element with given key and node pointer
El         *El_new(const char *key, PopNode *node) {
    El         *new = malloc(sizeof(*new));
    CHECKMEM(new);

    new->next = NULL;
    snprintf(new->key, sizeof(new->key), "%s", key);
    assert(0 == strncmp(new->key, key, sizeof(new->key)));
    new->node = node;
    return new;
}

/// Insert a new key/node pair into the linked list. Usage:
///
///     El *list=NULL;
///     int status;
///     list = El_insert(list, key, node, &status);
///     if(status != 0)
///         printf("Error: key is already in list\n");
///
/// @param self current element of list
/// @param[in] key a character string, the name of the current PopNode
/// @param[in] node pointer to PopNode object
/// @param[out] status pointer to int, which will be set to 0 on
/// success or 1 if node is already in list.
/// @return pointer to self or to a newly-allocated El.
El *El_insert(El *self, const char *key, PopNode *node, int *status) {
    El *new;
    if(self == NULL) {
        new = El_new(key, node);
        new->next = NULL;
        *status = 0; // success
        return new;
    }
    int cmp = strncmp(key, self->key, sizeof(new->key));
    if(cmp < 0) {
        new = El_new(key, node);
        new->next = self;
        *status = 0; // success
        return new;
    }else if(cmp > 0) {
        self->next = El_insert(self->next, key, node, status);
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

/// PopNodeTab constructor
PopNodeTab    *PopNodeTab_new(void) {
    PopNodeTab    *new = malloc(sizeof(*new));
    CHECKMEM(new);
    memset(new->tab, 0, sizeof(new->tab));
    return new;
}

/// PopNodeTab destructor
void PopNodeTab_free(PopNodeTab * self) {
    int i;
    for(i=0; i < PNT_DIM; ++i)
        El_free(self->tab[i]);
    free(self);
}

/// Get the PopNode object associated with key. If no
/// such object exists, return NULL.
PopNode *PopNodeTab_get(PopNodeTab * self, const char *key) {

    // Same as hash % PNT_DIM but faster. Requires
    // that PNT_DIM be a power of 2.
    unsigned h = strhash(key) & (PNT_DIM-1u);
    assert(h < PNT_DIM);

    assert(self);

    El *el;
    for(el = self->tab[h]; el; el = el->next) {
        int cmp = strncmp(key, el->key, sizeof(el->key));
        if(cmp == 0)
            return el->node;
        else if(cmp < 0)
            return NULL;
    }
    return NULL;
}

/// Insert a PopNode into the table.
/// @return 0 on success; 1 if PopNode is already in table.
int PopNodeTab_insert(PopNodeTab * self, const char *key,
                      PopNode *node) {

    // Same as hash % PNT_DIM but faster. Requires
    // that PNT_DIM be a power of 2.
    unsigned h = strhash(key) & (PNT_DIM-1u);
    assert(h < PNT_DIM);
    int status;

    assert(self);

    self->tab[h] = El_insert(self->tab[h], key, node, &status);
    return status;
}

/// Return the number of elements in the PopNodeTab.
unsigned long PopNodeTab_size(PopNodeTab * self) {
    unsigned    i;
    unsigned long size = 0;

    for(i = 0; i < PNT_DIM; ++i) {
        El *el;
        for(el = self->tab[i]; el; el = el->next)
            ++size;
    }
    return size;
}

/// Print a PopNodeTab
void PopNodeTab_print(PopNodeTab *self) {
    unsigned i;
    for(i=0; i < PNT_DIM; ++i) {
        printf("%2u:", i);
        El *el;
        for(el = self->tab[i]; el; el = el->next)
            printf(" [%s, %p]", el->key, el->node);
        putchar('\n');
    }
}

/// Return root of population tree. Run PopNodeTab_sanityCheck before
/// calling this function.
PopNode *PopNodeTab_root(PopNodeTab *self) {
    int i;
    for(i=0; i < PNT_DIM; ++i) {
        if(self->tab[i])
            return PopNode_root(self->tab[i]->node);
    }
    return NULL;
}

/// Check the sanity of each node and make sure there is only one root.
void PopNodeTab_sanityCheck(PopNodeTab *self, const char *file, int line) {
    int i;
    PopNode *node, *root = NULL;
    for(i=0; i < PNT_DIM; ++i) {
        El *el;
        for(el = self->tab[i]; el; el = el->next) {
            PopNode_sanityFromLeaf(el->node, file, line);
            node = PopNode_root(el->node);
            if(root == NULL)
                root = node;
            else if(root != node) {
                fprintf(stderr,
                        "%s:%d: Pop tree has multiple roots.\n",
                        file, line);
                exit(EXIT_FAILURE);
            }
        }
    }
}

