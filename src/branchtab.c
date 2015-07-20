#include "branchtab.h"
#include "misc.h"
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// BT_DIM must be a power of 2
#define BT_DIM 16u

// Make sure BT_DIM is a power of 2
#if (BT_DIM==0u || (BT_DIM & (BT_DIM-1u)))
# error BT_DIM must be a power of 2
#endif

// A single element
struct BTLink {
    struct BTLink  *next;
    unsigned        key;
    double          value;
};

struct BranchTab {
    BTLink         *tab[BT_DIM];
};

unsigned    typeIdHash(typeId_t key);
BTLink     *BTLink_add(typeId_t key, double value);
void        BTLink_free(BTLink * e);

BTLink         *BTLink_new(const char *key) {
    BTLink         *new = malloc(sizeof(*new));
    CHECKMEM(new);

    new->next = NULL;
    snprintf(new->key, sizeof(new->key), "%s", key);
    assert(0 == strncmp(new->key, key, sizeof(new->key)));
    new->value = NULL;
    return new;
}

void BTLink_free(BTLink * e) {
    if(e == NULL)
        return;
    BTLink_free(e->next);
    free(e);
}

void * BTLink_get(BTLink * self) {
    assert(self);
    return self->value;
}

void BTLink_set(BTLink * self, void *value) {
    assert(self);
    self->value = value;
}

BTLink *BTLink_find(BTLink * self, BTLink ** found, const char *key) {
    int comparison;
    if(self == NULL
       || 0 > (comparison = strncmp(key, self->key, sizeof(self->key)))) {
        *found = BTLink_new(key);
        (*found)->next = self;
        return *found;
    } else if(comparison > 0) {
        self->next = BTLink_find(self->next, found, key);
        return self;
    }
    assert(0 == strncmp(key, self->key, sizeof(self->key)));
    *found = self;
    return self;
}

void BTLink_printShallow(BTLink * self) {
    if(self == NULL)
        return;
    printf(" [%s, %p]", self->key, self->value);
}

void BTLink_print(BTLink * self) {
    if(self == NULL)
        return;
    BTLink_printShallow(self);
    BTLink_print(self->next);
}

BranchTab    *BranchTab_new(void) {
    BranchTab    *new = malloc(sizeof(*new));
    CHECKMEM(new);
    memset(new->tab, 0, sizeof(new->tab));
    return new;
}

void BranchTab_free(BranchTab * self) {
    free(self);
}

BTLink *BranchTab_get(BranchTab * self, const char *key) {

    unsigned h = strhash(key);
    assert(h < BT_DIM);

    assert(self);

    BTLink *el = NULL;
    self->tab[h] = BTLink_find(self->tab[h], &el, key);
    assert(el);

    return el;
}

/// Return the number of elements in the BranchTab.
unsigned long BranchTab_size(BranchTab * self) {
    unsigned    i;
    unsigned long size = 0;

    for(i = 0; i < BT_DIM; ++i) {
        BTLink *el;
        for(el = self->tab[i]; el; el = el->next)
            ++size;
    }
    return size;
}

void BranchTab_print(BranchTab *self) {
    unsigned i;
    for(i=0; i < BT_DIM; ++i) {
        printf("%2u:", i);
        BTLink_print(self->tab[i]);
        putchar('\n');
    }
}

