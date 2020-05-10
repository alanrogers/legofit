/**
 * @file strdblmap.c
 * @author Alan R. Rogers
 * @brief Associate character strings with doubles.
 *
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "strdblmap.h"
#include "misc.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

/// STRDBL_DIM must be a power of 2
#define STRDBL_DIM 64u

#if (STRDBL_DIM==0u || (STRDBL_DIM & (STRDBL_DIM-1u)))
#  error STRDBL_DIM must be a power of 2
#endif

#define MAXKEY 10

/// A single element of a linked list
typedef struct SDLink {
    struct SDLink *next;
    char        key[MAXKEY];
    double      value;
} SDLink;

/// Hash table
struct StrDblMap {
    SDLink     *tab[STRDBL_DIM];
};

SDLink     *SDLink_new(const char *key, double value, SDLink * next);
void        SDLink_free(SDLink * self);
SDLink     *SDLink_insert(SDLink * self, const char *key, double value);
double      SDLink_get(SDLink * self, const char *key);
int         SDLink_exists(SDLink * self, const char *key);
void        SDLink_print(const SDLink * self, FILE *fp);
unsigned    SDLink_size(SDLink *self);

/// SDLink constructor
SDLink     *SDLink_new(const char *key, double value, SDLink * next) {
    SDLink     *new = malloc(sizeof(*new));
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
unsigned SDLink_size(SDLink *self) {
    if(self == NULL)
        return 0u;
    return 1u + SDLink_size(self->next);
}

/// Free linked list of SDLink objects
void SDLink_free(SDLink * self) {
    if(self == NULL)
        return;
    SDLink_free(self->next);
    free(self);
}

/// Insert a new key-value pair. Set errno=EDOM if key already exists.
SDLink     *SDLink_insert(SDLink * self, const char *key, double value) {
    if(self == NULL)
        return SDLink_new(key, value, self);

    int         diff = strcmp(key, self->key);
    if(diff == 0) {
		errno = EDOM;
		return self;
    }else if(diff > 0) {
        self->next = SDLink_insert(self->next, key, value);
        return self;
    }else
        return SDLink_new(key, value, self);
}

/// Get index corresponding to key. If key is not in table, return -1 and
/// set errno = EDOM.
double SDLink_get(SDLink * self, const char *key) {
    if(self == NULL) {
		errno = EDOM;
		return -1.0;
    }

    int diff = strcmp(key, self->key);
    if(diff == 0)
        return self->value;
    else if(diff > 0)
        return SDLink_get(self->next, key);
    else{
        assert(diff < 0);
		errno = EDOM;
        return -1.0;
    }
}

/// Return 1 if key is in list, 0 otherwise
int SDLink_exists(SDLink * self, const char *key) {
    if(self == NULL)
		return 0;

    int diff = strcmp(key, self->key);
    if(diff == 0)
        return 1;
    else if(diff > 0)
        return SDLink_exists(self->next, key);
    else{
        assert(diff < 0);
        return 0;
    }
}

/// Print linked list of SDLink objects
void SDLink_print(const SDLink * self, FILE *fp) {
    if(self == NULL)
        return;
    fprintf(fp, " [%s, %g]", self->key, self->value);
    SDLink_print(self->next, fp);
}

/// StrDblMap constructor
StrDblMap     *StrDblMap_new(void) {
    StrDblMap     *new = malloc(sizeof(*new));
    CHECKMEM(new);
    memset(new, 0, sizeof(*new));
    return new;
}

/// StrDblMap destructor
void StrDblMap_free(StrDblMap * self) {
    int         i;
    for(i = 0; i < STRDBL_DIM; ++i)
        SDLink_free(self->tab[i]);
    free(self);
}

/// Insert a key-value pair into the hash table. Set errno=EDOM if
/// pair already exists.
void StrDblMap_insert(StrDblMap *self, const char *key, double value) {
    unsigned long h = strhash(key) & (STRDBL_DIM - 1ul);
    assert(h < STRDBL_DIM);
    assert(self);
    self->tab[h] = SDLink_insert(self->tab[h], key, value);
}

/// Return value corresponding to key. If key is not in table, return
/// -1 and set errno = EDOM.
double StrDblMap_get(StrDblMap * self, const char *key) {
    unsigned long h = strhash(key) & (STRDBL_DIM - 1ul);
    assert(h < STRDBL_DIM);
    assert(self);
    return SDLink_get(self->tab[h], key);
}

int StrDblMap_exists(StrDblMap * self, const char *key) {
    unsigned long h = strhash(key) & (STRDBL_DIM - 1ul);
    assert(h < STRDBL_DIM);
    assert(self);
    return SDLink_exists(self->tab[h], key);
}

/// Print a StrDblMap object
void StrDblMap_print(const StrDblMap * self, FILE *fp) {
    unsigned    i;
    for(i = 0; i < STRDBL_DIM; ++i) {
        fprintf(fp, "%2u:", i);
        SDLink_print(self->tab[i], fp);
        putc('\n', fp);
    }
}

/// Number of items stored in hash table.
unsigned StrDblMap_size(const StrDblMap *self) {
    unsigned i, n=0;

    for(i=0; i < STRDBL_DIM; ++i)
        n += SDLink_size(self->tab[i]);
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
            fprintf(stderr, "usage: xstrdbl [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    int         i;
    char        key[20];
    StrDblMap     *sd = StrDblMap_new();
    CHECKMEM(sd);
    assert(0 == StrDblMap_size(sd));

    for(i = 0; i < 100; ++i) {
        snprintf(key, sizeof key, "%d", i);
        StrDblMap_insert(sd, key, (double) i);
    }
    assert(100 == StrDblMap_size(sd));
	errno = 0;
	StrDblMap_insert(sd, "1", 1);
	assert(errno == EDOM);
    for(i = 0; i < 100; ++i) {
        snprintf(key, sizeof key, "%d", i);
        assert(((double) i) == StrDblMap_get(sd, key));
    }
	errno = 0;
    fflush(stdout);
	assert(-1 == StrDblMap_get(sd, "notthere"));
	assert(errno == EDOM);
    assert(0 == StrDblMap_exists(sd, "notthere"));
    assert(1 == StrDblMap_exists(sd, "1"));

    if(verbose)
        StrDblMap_print(sd, stdout);

    StrDblMap_free(sd);

    unitTstResult("StrDblMap", "OK");
}
#endif
