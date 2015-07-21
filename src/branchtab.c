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
typedef struct BTLink {
    struct BTLink  *next;
    tipId_t         key;
    double          value;
} BTLink;

struct BranchTab {
    BTLink         *tab[BT_DIM];
};

BTLink     *BTLink_new(tipId_t key, double value);
BTLink     *BTLink_add(BTLink * self, tipId_t key, double value);
double      BTLink_get(BTLink * self, tipId_t key);
void        BTLink_free(BTLink * self);
void        BTLink_printShallow(BTLink * self);
void        BTLink_print(BTLink * self);

#if TIPID_SIZE==32
uint32_t    tipIdHash(uint32_t key);
#elif TIPID_SIZE==64
uint32_t    tipIdHash(uint64_t key);
#endif

#if TIPID_SIZE==32
// These hash functions are from Thomas Wang's 1997 article:
// https://gist.github.com/badboy/6267743
//
// This one hashes a 32-bit integer 
uint32_t tipIdHash( uint32_t key) {
   key = (key+0x7ed55d16) + (key<<12);
   key = (key^0xc761c23c) ^ (key>>19);
   key = (key+0x165667b1) + (key<<5);
   key = (key+0xd3a2646c) ^ (key<<9);
   key = (key+0xfd7046c5) + (key<<3);
   key = (key^0xb55a4f09) ^ (key>>16);
   return key & (BT_DIM-1);
}
#elif TIPID_SIZE==64
// This one hashes a 64-bit integer but returns 32 bits.
uint32_t tipIdHash(uint64_t key) {
  key = (~key) + (key << 18); // key = (key << 18) - key - 1;
  key = key ^ (key >> 31);
  key = key * 21; // key = (key + (key << 2)) + (key << 4);
  key = key ^ (key >> 11);
  key = key + (key << 6);
  key = key ^ (key >> 22);
  return (uint32_t) (key  & (BT_DIM-1));
}
#else
#error "Can't compile tipIdHash function. See branchtab.c"
#endif

BTLink         *BTLink_new(tipId_t key, double value) {
    BTLink         *new = malloc(sizeof(*new));
    CHECKMEM(new);

    new->next = NULL;
    new->key = key;
    new->value = value;
    return new;
}

void BTLink_free(BTLink * self) {
    if(self == NULL)
        return;
    BTLink_free(self->next);
    free(self);
}

BTLink *BTLink_add(BTLink * self, tipId_t key, double value) {
    if(self == NULL || key < self->key) {
        BTLink *new = BTLink_new(key, value);
        new->next = self;
        return new;
    } else if(key > self->key) {
        self->next = BTLink_add(self->next, key, value);
        return self;
    }
    assert(key == self->key);
    self->value += value;
    return self;
}

/// Return value corresponding to key, or nan if no value is found.
double BTLink_get(BTLink * self, tipId_t key) {
    if(self == NULL || key < self->key)
        return nan("");
    else if(key > self->key)
        return BTLink_get(self->next, key);
    assert(key == self->key);
    return self->value;
}

void BTLink_printShallow(BTLink * self) {
    if(self == NULL)
        return;
    printf(" [%lu, %lf]", (unsigned long) self->key, self->value);
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
    memset(new, 0, sizeof(*new));
    return new;
}

void BranchTab_free(BranchTab * self) {
    int i;
    for(i=0; i < BT_DIM; ++i)
        BTLink_free(self->tab[i]);
    free(self);
}

/// Return value corresponding to key, or nan if no value is found.
double BranchTab_get(BranchTab * self, tipId_t key) {
    unsigned h = tipIdHash(key);
    assert(h < BT_DIM);
    assert(self);
    return BTLink_get(self->tab[h], key);
}

/// Add a value to table. If key already exists, new value is added to
/// old one.  
void BranchTab_add(BranchTab * self, tipId_t key, double value) {
    unsigned h = tipIdHash(key);
    assert(h < BT_DIM);
    assert(self);
    self->tab[h] = BTLink_add(self->tab[h], key, value);
}

/// Return the number of elements in the BranchTab.
unsigned BranchTab_size(BranchTab * self) {
    unsigned    i;
    unsigned    size = 0;

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

/// Add each entry in table rhs to table lhs
void BranchTab_plusEquals(BranchTab *lhs, BranchTab *rhs) {
    int i;
    for(i=0; i<BT_DIM; ++i) {
        BTLink *link;
        for(link=rhs->tab[i]; link!=NULL; link=link->next)
            BranchTab_add(lhs, link->key, link->value);
    }
}

void BranchTab_toArrays(BranchTab *self, unsigned n, tipId_t key[n], double value[n]) {
    int i, j=0;
    for(i=0; i<BT_DIM; ++i) {
        BTLink *link;
        for(link=self->tab[i]; link!=NULL; link=link->next) {
            if(j >= n)
                eprintf("%s:%s:%d: buffer overflow\n", __FILE__,__func__,__LINE__);
            key[j] = link->key;
            value[j] = link->value;
            ++j;
        }
    }
}

#ifdef TEST

#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {

    int verbose=0;

    if(argc > 1) {
        if(argc!=2 || 0!=strcmp(argv[1], "-v")) {
            fprintf(stderr,"usage: xbranchtab [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    BranchTab *bt = BranchTab_new();
    CHECKMEM(bt);
    assert(0 == BranchTab_size(bt));

    tipId_t key[25];
    double  val[25];

    int i;
    for(i=0; i < 25; ++i) {
        key[i] = i+1;
        val[i] = (double) key[i];
        BranchTab_add(bt, key[i], val[i]);
        assert(i+1 == BranchTab_size(bt));
    }
    for(i=0; i < 25; ++i) {
        key[i] = i+1;
        val[i] = (double) key[i];
        BranchTab_add(bt, key[i], val[i]);
    }

    for(i=0; i < 25; ++i) {
        assert(2*val[i] == BranchTab_get(bt, key[i]));
    }

    if(verbose)
        BranchTab_print(bt);
    BranchTab_free(bt);
}
#endif
