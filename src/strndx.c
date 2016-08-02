/**
 * @file strndx.c
 * @author Alan R. Rogers
 * @brief Associate character strings with a 0-based index.
 *
 * @copyright Copyright (c) 2016, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "strndx.h"
#include "misc.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

// STRNDX_DIM must be a power of 2
#define STRNDX_DIM 64u

#if (STRNDX_DIM==0u || (STRNDX_DIM & (STRNDX_DIM-1u)))
#  error STRNDX_DIM must be a power of 2
#endif

#define MAXKEY 10

// A single element
typedef struct SNLink {
    struct SNLink *next;
    char        key[MAXKEY];
    int         ndx;
} SNLink;

struct StrNdx {
    int         nitems;
    SNLink     *tab[STRNDX_DIM];
};

SNLink     *SNLink_new(char *key, int ndx, SNLink * next);
SNLink     *SNLink_get(SNLink * self, char *key, int *ndx, int *nitems);
void        SNLink_free(SNLink * self);
void        SNLink_print(const SNLink * self);

SNLink     *SNLink_new(char *key, int ndx, SNLink * next) {
    SNLink     *new = malloc(sizeof(*new));
    CHECKMEM(new);

    new->next = next;
    int         status = snprintf(new->key, sizeof new->key, "%s", key);
    if(status >= MAXKEY) {
        fprintf(stderr, "%s:%s:%d: Buffer overflow. MAXKEY=%d, key=%s\n",
                __FILE__, __func__, __LINE__, MAXKEY, key);
        free(new);
        exit(EXIT_FAILURE);
    }

    new->ndx = ndx;
    return new;
}

void SNLink_free(SNLink * self) {
    if(self == NULL)
        return;
    SNLink_free(self->next);
    free(self);
}

/// Get index corresponding to key. If key has never been seen before,
/// a new key-index item is created and placed into the hash
/// table. If the key already exists, its index is placed in *ndx.
SNLink     *SNLink_get(SNLink * self, char *key, int *ndx, int *nitems) {
    if(self == NULL) {
        *ndx = *nitems;
        ++*nitems;
        return SNLink_new(key, *ndx, self);
    }

    int         diff = strcmp(key, self->key);
    if(diff == 0) {
        *ndx = self->ndx;
        return self;
    }else if(diff > 0) {
        self->next = SNLink_get(self->next, key, ndx, nitems);
        return self;
    }else{
        assert(diff < 0);
        *ndx = *nitems;
        ++*nitems;
        return SNLink_new(key, *ndx, self);
    }
}

void SNLink_print(const SNLink * self) {
    if(self == NULL)
        return;
    printf(" [%s, %d]", self->key, self->ndx);
    SNLink_print(self->next);
}

StrNdx     *StrNdx_new(void) {
    StrNdx     *new = malloc(sizeof(*new));
    CHECKMEM(new);
    memset(new, 0, sizeof(*new));
    return new;
}

void StrNdx_free(StrNdx * self) {
    int         i;
    for(i = 0; i < STRNDX_DIM; ++i)
        SNLink_free(self->tab[i]);
    free(self);
}

/// Return ndx corresponding to key; abort on failure
int StrNdx_getNdx(StrNdx * self, char *key) {
    int ndx;
    unsigned    h = strhash(key) & (STRNDX_DIM - 1u);
    assert(h < STRNDX_DIM);
    assert(self);
    self->tab[h] = SNLink_get(self->tab[h], key, &ndx, &self->nitems);
    return ndx;
}

/// Return the number of elements in the StrNdx.
int StrNdx_size(StrNdx * self) {
    return self->nitems;
}

void StrNdx_print(const StrNdx * self) {
    unsigned    i;
    for(i = 0; i < STRNDX_DIM; ++i) {
        printf("%2u:", i);
        SNLink_print(self->tab[i]);
        putchar('\n');
    }
}

#ifdef TEST

#  include <string.h>
#  include <assert.h>
#  include <unistd.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

int main(int argc, char **argv) {
    int         verbose = 0;
    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xstrndx [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    StrNdx     *sn = StrNdx_new();
    CHECKMEM(sn);
    assert(0 == StrNdx_size(sn));

    char        key[20];

    int         i, j;
    for(i = 0; i < 10; ++i) {
        snprintf(key, sizeof key, "%d", i);
        j = StrNdx_getNdx(sn, key);
        assert(i==j);
    }

    assert(10 == StrNdx_size(sn));

    for(i = 0; i < 25; ++i) {
        snprintf(key, sizeof key, "%d", i);
        assert(i==StrNdx_getNdx(sn, key));
    }

    if(verbose)
        StrNdx_print(sn);
    assert(25 == StrNdx_size(sn));

    StrNdx_free(sn);

    unitTstResult("StrNdx", "OK");
}
#endif
