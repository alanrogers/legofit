/**
 * @file strint.c
 * @author Alan R. Rogers
 * @brief Associate character strings with integers.
 *
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "strint.h"
#include "misc.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

/// STRINT_DIM must be a power of 2
#define STRINT_DIM 64u

#if (STRINT_DIM==0u || (STRINT_DIM & (STRINT_DIM-1u)))
#  error STRINT_DIM must be a power of 2
#endif

#define MAXKEY 10

/// A single element of a linked list
typedef struct SILink {
    struct SILink *next;
    char        key[MAXKEY];
    int         value;
} SILink;

/// Hash table
struct StrInt {
    SILink     *tab[STRINT_DIM];
};

SILink     *SILink_new(const char *key, int value, SILink * next);
void        SILink_free(SILink * self);
SILink     *SILink_insert(SILink * self, const char *key, int value);
int         SILink_get(SILink * self, const char *key);
void        SILink_print(const SILink * self, FILE *fp);
unsigned    SILink_size(SILink *self);

/// SILink constructor
SILink     *SILink_new(const char *key, int value, SILink * next) {
    SILink     *new = malloc(sizeof(*new));
    CHECKMEM(new);

    new->next = next;
    new->value = value;
    int         status = snprintf(new->key, sizeof new->key, "%s", key);
    if(status >= MAXKEY) {
        fprintf(stderr, "%s:%s:%d: Buffer overflow. MAXKEY=%d, key=%s\n",
                __FILE__, __func__, __LINE__, MAXKEY, key);
        free(new);
        exit(EXIT_FAILURE);
    }

    return new;
}

/// Return number of links in list
unsigned SILink_size(SILink *self) {
    if(self == NULL)
        return 0u;
    return 1u + SILink_size(self->next);
}

/// Free linked list of SILink objects
void SILink_free(SILink * self) {
    if(self == NULL)
        return;
    SILink_free(self->next);
    free(self);
}

/// Insert a new key-value pair. Set errno=EDOM if key already exists.
SILink     *SILink_insert(SILink * self, const char *key, int value) {
    if(self == NULL)
        return SILink_new(key, value, self);

    int         diff = strcmp(key, self->key);
    if(diff == 0) {
		errno = EDOM;
		return self;
    }else if(diff > 0) {
        self->next = SILink_insert(self->next, key, value);
        return self;
    }else
        return SILink_new(key, value, self);
}

/// Get index corresponding to key. If key is not in table, return -1 and
/// set errno = EDOM.
int SILink_get(SILink * self, const char *key) {
    if(self == NULL) {
		errno = EDOM;
		return -1;
    }

    int diff = strcmp(key, self->key);
    if(diff == 0)
        return self->value;
    else if(diff > 0)
        return SILink_get(self->next, key);
    else{
        assert(diff < 0);
		errno = EDOM;
        return -1;
    }
}

/// Print linked list of SILink objects
void SILink_print(const SILink * self, FILE *fp) {
    if(self == NULL)
        return;
    fprintf(fp, " [%s, %d]", self->key, self->value);
    SILink_print(self->next, fp);
}

/// StrInt constructor
StrInt     *StrInt_new(void) {
    StrInt     *new = malloc(sizeof(*new));
    CHECKMEM(new);
    memset(new, 0, sizeof(*new));
    return new;
}

/// StrInt destructor
void StrInt_free(StrInt * self) {
    int         i;
    for(i = 0; i < STRINT_DIM; ++i)
        SILink_free(self->tab[i]);
    free(self);
}

/// Insert a key-value pair into the hash table. Set errno=EDOM if
/// pair already exists.
void StrInt_insert(StrInt *self, const char *key, int value) {
    unsigned    h = strhash(key) & (STRINT_DIM - 1u);
    assert(h < STRINT_DIM);
    assert(self);
    self->tab[h] = SILink_insert(self->tab[h], key, value);
}

/// Return value corresponding to key. If key is not in table, return
/// -1 and set errno = EDOM.
int StrInt_get(StrInt * self, const char *key) {
    unsigned    h = strhash(key) & (STRINT_DIM - 1u);
    assert(h < STRINT_DIM);
    assert(self);
    return SILink_get(self->tab[h], key);
}

/// Print a StrInt object
void StrInt_print(const StrInt * self, FILE *fp) {
    unsigned    i;
    for(i = 0; i < STRINT_DIM; ++i) {
        fprintf(fp, "%2u:", i);
        SILink_print(self->tab[i], fp);
        putc('\n', fp);
    }
}

/// Number of items stored in hash table.
unsigned StrInt_size(const StrInt *self) {
    unsigned i, n=0;

    for(i=0; i < STRINT_DIM; ++i)
        n += SILink_size(self->tab[i]);
    return n;
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
            fprintf(stderr, "usage: xstrint [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    int         i;
    char        key[20];
    StrInt     *si = StrInt_new();
    CHECKMEM(si);
    assert(0 == StrInt_size(si));

    for(i = 0; i < 100; ++i) {
        snprintf(key, sizeof key, "%d", i);
        StrInt_insert(si, key, i);
    }
    assert(100 == StrInt_size(si));
	errno = 0;
	StrInt_insert(si, "1", 1);
	assert(errno == EDOM);
    for(i = 0; i < 100; ++i) {
        snprintf(key, sizeof key, "%d", i);
        assert(i == StrInt_get(si, key));
    }
	errno = 0;
    fflush(stdout);
	assert(-1 == StrInt_get(si, "notthere"));
	assert(errno == EDOM);

    if(verbose)
        StrInt_print(si, stdout);

    StrInt_free(si);

    unitTstResult("StrInt", "OK");
}
#endif
