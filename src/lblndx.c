/**
 * @file lblndx.c
 * @author Alan R. Rogers
 * @brief An index of sample labels.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "error.h"
#include "lblndx.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

int         comparePtrs(const void *void_x, const void *void_y);
tipId_t     LblNdx_getTipId_1(const LblNdx *self, const char *lbl);

/// Set everything to zero.
void LblNdx_init(LblNdx * self) {
    memset(self, 0, sizeof(*self));
    LblNdx_sanityCheck(self, __FILE__, __LINE__);
}

/// Add samples for a single population. Should be called once for
/// each sampled population.
void LblNdx_addSamples(LblNdx * self, unsigned nsamples, const char *lbl) {
    unsigned    i;
    if(self->n + nsamples >= MAXSAMP)
        eprintf("%s:%s:%d: too many samples\n", __FILE__, __func__, __LINE__);
    for(i = 0; i < nsamples; ++i) {
        if(nsamples == 1)
            snprintf(self->lbl[self->n], POPNAMESIZE, "%s", lbl);
        else
            snprintf(self->lbl[self->n], POPNAMESIZE, "%s.%u", lbl, i);
        self->n += 1;
    }
    LblNdx_sanityCheck(self, __FILE__, __LINE__);
}

/// Return the label associated with index i.
const char *LblNdx_lbl(const LblNdx * self, unsigned i) {
    assert(i < self->n);
    return self->lbl[i];
}

/// Return number of samples, which equals number of labels.
unsigned LblNdx_size(const LblNdx * self) {
    return self->n;
}

void        LblNdx_sanityCheck(const LblNdx *self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self, file, line);
    REQUIRE(self->n < MAXSAMP, file, line);

    int i;
    for(i=0; i < self->n; ++i) {
        REQUIRE(NULL != self->lbl[i], file, line);
        REQUIRE(legalName(self->lbl[i]), file, line);
    }
#endif
}

/// Return 1 if the two arguments are equal; 0 otherwise.
int         LblNdx_equals(const LblNdx *lhs, const LblNdx *rhs) {
    if(lhs == rhs)
        return 1;
    if(lhs==NULL || rhs==NULL)
        return 0;
    if(lhs->n != rhs->n)
        return 0;

    int i;
    for(i=0; i < lhs->n; ++i)
        if(0 != strcmp(lhs->lbl[i], rhs->lbl[i]))
            return 0;
    return 1;
}

// Return tipId_t value corresponding to a label. The label should
// refer to a single population. In other words, it should not contain
// the ":" character. For the i'th label, this value equals the i'th
// power of 2.
tipId_t     LblNdx_getTipId_1(const LblNdx *self, const char *lbl) {
    unsigned i;
    const tipId_t unity = 1;

    for(i=0; i < self->n; ++i) {
        if(0 == strcmp(lbl, self->lbl[i]))
           break;
    }
    if(i == self->n)
        return 0;

    return unity << i;
}

/// Reverse lookup. Return tipId_t value corresponding to
/// label, which may be composite. The parts of a composite
/// label are delimited by the ':' character. If label (or any of its
/// parts) is not present in LblNdx, return 0. 
tipId_t LblNdx_getTipId(const LblNdx *self, const char *lbl) {
    char buff[200];
    const char *start=lbl, *end;
    int len;
    tipId_t rval = 0, id;
    while(1) {
        end = strchr(start, ':');

        // Not a composite label
        if(end == NULL) {
            id = LblNdx_getTipId_1(self, start);
            if(id == 0)
                return 0;
            rval |= id;
            break; // loop exit
        }

        // composite label: copy first component into
        // buff, get id of that component, and "or" the
        // component id into the return value.
        len = end - start;
        int status = strnncopy(sizeof(buff), buff, len, start);
        switch(status) {
        case 0:
            break;
        case BUFFER_OVERFLOW:
            fprintf(stderr,"%s:%d: buffer overflow\n", __FILE__,__LINE__);
            exit(EXIT_FAILURE);
        default:
            fprintf(stderr,"%s:%d: unknown error\n", __FILE__,__LINE__);
            exit(EXIT_FAILURE);
        }
        id = LblNdx_getTipId_1(self, buff);
        if(id == 0)
            return 0;
        rval |= id;
        start = end+1;
    }
    return rval;
}

/// Print a LblNdx
void LblNdx_print(const LblNdx *self, FILE *fp) {
    unsigned i;
    for(i=0; i < self->n; ++i)
        fprintf(fp,"%3u %5s\n", i, self->lbl[i]);
}

/// Generate a label for site pattern tid. Label goes into
/// buff. Function returns a pointer to buff;
char       *patLbl(size_t n, char buff[n], tipId_t tid, const LblNdx * lblndx) {
    const int   maxbits = 40;
    int         bit[maxbits];
    int         i, nbits;
    char        lbl[100];

    if(TIPID_SIZE - nlz(tid) > lblndx->n) {
        fprintf(stderr,"%s:%d: LblNdx only allows for %d bits;"
                " tid has %d.\n",
                __FILE__,__LINE__, lblndx->n,
                TIPID_SIZE - nlz(tid));
        exit(EXIT_FAILURE);
    }
    nbits = getBits(tid, maxbits, bit);
    buff[0] = '\0';
    for(i = 0; i < nbits; ++i) {
        snprintf(lbl, sizeof(lbl), "%s",
                 LblNdx_lbl(lblndx, (unsigned) bit[i]));
        if(strlen(buff) + strlen(lbl) >= n)
            eprintf("%s:%s:%d: buffer overflow\n", __FILE__, __func__,
                    __LINE__);
        strcat(buff, lbl);
        if(i + 1 < nbits && 1 + strlen(buff) < n)
            strcat(buff, ":");
    }
    return buff;
}

/// Compare pointers to pointers to two tipId_t values.
///
/// @param void_x,void_y pointers to pointers to tipId_t values
/// @returns <0, 0, or >0 depending on whether the first arg is <,
/// ==, or > the second.
int comparePtrs(const void *void_x, const void *void_y) {
    tipId_t * const * x = (tipId_t * const *) void_x;
    tipId_t * const * y = (tipId_t * const *) void_y;

    // Major sort is on the number of samples
    // represented in the site pattern. Patterns with
    // fewer samples come first.
    int diff1bits = num1bits(**x) - num1bits(**y);
    if(diff1bits)
        return diff1bits;

    // Reverse order of bits so that low-order bit
    // is most significant. This ensures that the
    // sort order of samples corresponds to the
    // order in which they were listed in the input
    // data.
    unsigned rx = reverseBits(**x);
    unsigned ry = reverseBits(**y);

    return ry>rx ? 1 : ry<rx ? -1 : 0;
}

/// Compare pointers to two tipId_t values.
///
/// @param void_x,void_y pointers to tipId_t values
/// @returns <0, 0, or >0 depending on whether the first arg is <,
/// ==, or > the second.
int compare_tipId(const void *void_x, const void *void_y) {
    tipId_t const * x = (tipId_t const *) void_x;
    tipId_t const * y = (tipId_t const *) void_y;

    // Major sort is on the number of samples
    // represented in the site pattern. Patterns with
    // fewer samples come first.
    int diff1bits = num1bits(*x) - num1bits(*y);
    if(diff1bits)
        return diff1bits;

    // Reverse order of bits so that low-order bit
    // is most significant. This ensures that the
    // sort order of samples corresponds to the
    // order in which they were listed in the input
    // data.
    unsigned rx = reverseBits(*x);
    unsigned ry = reverseBits(*y);

    return ry>rx ? 1 : ry<rx ? -1 : 0;
}

/// On entry, pat is an array of n tipId_t values. On return,
/// ord[0] is the index of the first value, ord[1] is that of the 2nd,
/// and so on.
void orderpat(int n, unsigned ord[n], tipId_t pat[n]) {
    tipId_t  *ptr[n];
    int i;
    for(i=0; i < n; ++i)
        ptr[i] = pat+i;
    qsort(ptr, (size_t) n, sizeof(ptr[0]), comparePtrs);
    for(i=0; i<n; ++i)
        ord[i] = ptr[i]-pat;
}

// Make map, an array whose i'th entry is the index in the new LblNdx
// of the i'th entry in the old LblNdx. Returns the index of the bit
// into which all "on" bits of "collapse" are mapped.
static int make_map(size_t n, int map[n], tipId_t collapse) {
    int i, min = n, shift=0;
    tipId_t bit = 1u;
    for(i=0; i < n; ++i, bit <<= 1) {
        if( collapse & bit ) {
            if(min == n) {
                min = i;
            }else
                ++shift;
            map[i] = min;
        }else
            map[i] = i - shift;
    }
    return min;
}

/**
 * Reduce dimension of LblNdx object by collapsing several entries
 * into a single entry. The "on" bits of "collapse" indicate the
 * entries to be collapsed.  they are collapsed into the entry
 * indicated by the bit in the lowest position.  "lbl" is the label
 * assigned to this entry. The dimension of the rewritten LblNdx
 * object is reduced by one less than the number of "on" bits in
 * "collapse".  The function returns 0 on success, EDOM if the number
 * of "on" bits in "collapse" exceeds self->n, and BUFFER_OVERFLOW if
 * any of the write operations would overflow.
 */
int LblNdx_collapse(LblNdx *self, tipId_t collapse, const char *lbl) {
    // Most significant set bit in collapse cannot be at position
    // greater than self->n.
    if(self->n < TIPID_SIZE - nlz(collapse))
        return EDOM;

    // Make map, an array whose i'th entry is the index in the new id
    // of the i'th bit in the old id.
    int map[self->n];
    int min = make_map(self->n, map, collapse);

    // Rewrite LblNdx object
    for(int i=0; i < self->n; ++i) {
        int j = map[i];
        if(j == min) {
            if(strlen(lbl) > POPNAMESIZE-1)
                return BUFFER_OVERFLOW;
            fprintf(stderr,"%s:%d:i=%d j=%d\n",
                    __FILE__,__LINE__,i,j);
            strcpy(self->lbl[j], lbl);
            fprintf(stderr,"%s:%d\n",__FILE__,__LINE__);
        }else {
            fprintf(stderr,"%s:%d:i=%d j=%d\n",
                    __FILE__,__LINE__,i,j);
            strcpy(self->lbl[j], self->lbl[i]);
        }
    }
    fprintf(stderr,"%s:%d\n",__FILE__,__LINE__);
    self->n -= num1bits(collapse) - 1;
    return 0;
}

#ifdef TEST

#  include <string.h>
#  include <assert.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

int main(int argc, char **argv) {

    int         verbose = 0;
    unsigned    i;

    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xlblndx [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    LblNdx     lndx = {.n = 3 };
    assert(lndx.n == 3);
    assert(LblNdx_size(&lndx) == 3);

    LblNdx_init(&lndx);
    assert(LblNdx_size(&lndx) == 0);

    LblNdx_addSamples(&lndx, 1, "A");
    LblNdx_addSamples(&lndx, 2, "B");

    assert(3 == LblNdx_size(&lndx));
    assert(0 == strcmp("A", LblNdx_lbl(&lndx, 0)));
    assert(0 == strcmp("B.0", LblNdx_lbl(&lndx, 1)));
    assert(0 == strcmp("B.1", LblNdx_lbl(&lndx, 2)));

    assert(01u == LblNdx_getTipId(&lndx, "A"));
    assert(02u == LblNdx_getTipId(&lndx, "B.0"));
    assert(04u == LblNdx_getTipId(&lndx, "B.1"));
    assert(03u == LblNdx_getTipId(&lndx, "A:B.0"));
    assert(05u == LblNdx_getTipId(&lndx, "A:B.1"));
    assert(06u == LblNdx_getTipId(&lndx, "B.0:B.1"));
    assert(07u == LblNdx_getTipId(&lndx, "A:B.0:B.1"));
    assert(07u == LblNdx_getTipId(&lndx, "B.0:A:B.1"));
    assert(0u == LblNdx_getTipId(&lndx, "x:A:B.1"));

    for(i=0; i < LblNdx_size(&lndx); ++i) {
        assert((1u << i) == LblNdx_getTipId(&lndx, LblNdx_lbl(&lndx, i)));
    }

    if(verbose)
        LblNdx_print(&lndx, stdout);

    LblNdx lndx2 = {.n = 3 };
    LblNdx_init(&lndx2);
    LblNdx_addSamples(&lndx2, 1, "A");
    LblNdx_addSamples(&lndx2, 2, "B");

    assert(LblNdx_equals(&lndx, &lndx2));

    // test LblNdx_collapse
    assert(3u == LblNdx_size(&lndx));
    tipId_t collapse = 06u; // collapse bits 2 and 3
    int status = LblNdx_collapse(&lndx, collapse, "D");
    switch(status) {
    case 0:
        break;
    case EDOM:
        fprintf(stderr,"%s:%d: bad input to LblNdx_collapse\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    case BUFFER_OVERFLOW:
        fprintf(stderr,"%s:%d: buffer overflow in LblNdx_collapse\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    default:
        fprintf(stderr,"%s:%d: unknown error in LblNdx_collapse\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    if(verbose)
        LblNdx_print(&lndx, stdout);
    assert(2 == LblNdx_size(&lndx));
    assert(0 == strcmp("A", LblNdx_lbl(&lndx, 0)));
    assert(0 == strcmp("D", LblNdx_lbl(&lndx, 1)));
    assert(1 == LblNdx_getTipId(&lndx, "A"));
    assert(2 == LblNdx_getTipId(&lndx, "D"));
    assert(3 == LblNdx_getTipId(&lndx, "A:D"));
           
    assert(3u == LblNdx_size(&lndx2));
    collapse = 011u; // illegal: bit too large
    status = LblNdx_collapse(&lndx, collapse, "D");
    assert(status == EDOM);

    unitTstResult("LblNdx_collapse", "OK");

    unitTstResult("LblNdx", "OK");

    return 0;
}
#endif
