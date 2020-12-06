/**
 * @file branchtab.c
 * @author Alan R. Rogers
 * @brief Hash table associating key (an unsigned int encoding a site
 * pattern) and value (a double representing the length of the
 * ascending branch)
 *
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "branchtab.h"
#include "misc.h"
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

struct BranchTab {
#ifndef NDEBUG
    int frozen;   // nonzero => no further changes allowed
#endif
    unsigned dim;
    unsigned nsamples;
    long double *tab; // zeroth entry is not used.
};

/// Create a new BranchTab with room for site patterns associated with
/// data in which the number of samples is "nsamples".
BranchTab    *BranchTab_new(unsigned nsamples) {
    if(nsamples < 2) {
        fprintf(stderr,"%s:%s:%d: nsamples=%u; must be > 1.\n",
                __FILE__,__func__,__LINE__, nsamples);
        exit(EXIT_FAILURE);
    }
    BranchTab    *new = malloc(sizeof(*new));
    CHECKMEM(new);

#ifndef NDEBUG
    new->frozen = 0;
#endif

    new->nsamples = nsamples;

    // Leave room for the site pattern in which all samples carry the
    // derived allele. This isn't ordinarily used, but can be produced
    // by BranchTab_collapse.
    new->dim = (1u << nsamples);
    
    new->tab = malloc(new->dim * sizeof(new->tab[0]));
    CHECKMEM(new->tab);

    memset(new->tab, 0, new->dim * sizeof(new->tab[0]));
    
    return new;
}

BranchTab    *BranchTab_dup(const BranchTab *old) {
    BranchTab *new = BranchTab_new(old->nsamples);
    CHECKMEM(new);

#ifndef NDEBUG
    new->frozen = old->frozen;
#endif

    new->nsamples = old->nsamples;
    new->dim = old->dim;
    new->tab = malloc(new->dim * sizeof(new->tab[0]));
    CHECKMEM(new->tab);

    memcpy(new->tab, old->tab, new->dim * sizeof(new->tab[0]));

    return new;
}

/// Return 1 if two BranchTab objects are equal; 0 otherwise.
int BranchTab_equals(const BranchTab *lhs, const BranchTab *rhs) {
#ifndef NDEBUG
    if(lhs->frozen != rhs->frozen)
        return 0;
#endif

    if(lhs->dim != rhs->dim)
        return 0;

    if(lhs->nsamples != rhs->nsamples)
        return 0;

    if(memcmp(lhs->tab, rhs->tab, lhs->dim * sizeof(lhs->tab[0])))
        return 0;
    return 1;
}

/// Destructor for BranchTab.
void BranchTab_free(BranchTab * self) {
    free(self->tab);
    free(self);
}

/// Return 1 if BranchTab includes singleton site patterns
int BranchTab_hasSingletons(BranchTab * self) {
    for(int i=1; i < self->dim; i <<= 1) {
        if(self->tab[i] != 0.0)
            return 1;
    }
    return 0;
}

/// Return value corresponding to key.
long double BranchTab_get(BranchTab * self, tipId_t key) {
    assert(key < self->dim);
    assert(key > 0);
    return self->tab[key];
}

/// Add a value to table. If key already exists, new value is added to
/// old one.
void BranchTab_add(BranchTab * self, tipId_t key, long double value) {
    assert(!self->frozen);
    assert(key > 0);
#ifndef NDEBUG    
    if(key >= self->dim) {
        fprintf(stderr,"%s:%d: key=o%o >= dim=o%o; nsamples=%u\n",
                __FILE__,__LINE__,key,self->dim,self->nsamples);
        exit(EXIT_FAILURE);
    }
#endif    
    assert(key < self->dim);
    self->tab[key] += value;
}

/// Return the number of elements in the BranchTab.
unsigned BranchTab_size(BranchTab * self) {
    return self->dim-1; // exclude self->tab[0]
}

void BranchTab_sanityCheck(BranchTab *self, const char *file, int line) {
#ifndef NDEBUG
    REQUIRE(self->nsamples > 1, file, line);
    REQUIRE(self->dim == (1u << self->nsamples), file, line);
    REQUIRE(self->tab != NULL, file, line);
    REQUIRE(self->tab[0] == 0, file, line);
    for(unsigned i=0; i < self->dim; ++i)
        REQUIRE( nlz(self->tab[i]) >= 8*sizeof(tipId_t) - self->nsamples,
                 file, line);
#endif    
}

/// Divide all values by denom. Return 0 on success, or 1 on failure.
int BranchTab_divideBy(BranchTab *self, long double denom) {
#ifndef NDEBUG
    assert(!self->frozen);
    self->frozen = 1;  // you can only call this function once
#endif

    // divide by denom
    for(unsigned i = 1; i < self->dim; ++i)
        self->tab[i] /= denom;

    return 0;
}

/// Print a BranchTab to standard output.
void BranchTab_print(const BranchTab *self, FILE *fp) {
    for(unsigned i=1; i < self->dim; ++i)
        fprintf(fp, "%u -> %Lg\n", i, self->tab[i]);
}

/// Add each entry in table rhs to table lhs
void BranchTab_plusEquals(BranchTab *lhs, BranchTab *rhs) {
    assert(!lhs->frozen);
    assert(lhs->dim == rhs->dim);
    for(unsigned i=1; i < lhs->dim; ++i)
        lhs->tab[i] += rhs->tab[i];
}

/// Subtract each entry in table rhs from table lhs
void BranchTab_minusEquals(BranchTab *lhs, BranchTab *rhs) {
    assert(!lhs->frozen);
    assert(lhs->dim == rhs->dim);
    for(unsigned i=1; i < lhs->dim; ++i)
        lhs->tab[i] -= rhs->tab[i];
}

/// Fill arrays key and value with values in BranchTab.  On return,
/// key[i] is the id of the i'th site pattern, and value[i] is the
/// total branch length associated with that site pattern.
void BranchTab_toArrays(BranchTab *self, unsigned n, tipId_t key[n],
                        long double value[n]) {
    assert(n >= self->dim-1);
    unsigned i;
    for(i=1; i < self->dim; ++i) {
        key[i-1] = i;
        value[i-1] = self->tab[i];
    }
}

/// Map two or more populations into a single population.
BranchTab *BranchTab_collapse(BranchTab *old, tipId_t collapse) {
    int newsamp = old->nsamples - num1bits(collapse) + 1;

    // Make map, an array whose i'th entry is an unsigned integer
    // with one bit on and the rest off. The on bit indicates the
    // position in the new id of the i'th bit in the old id.
    tipId_t map[old->nsamples];
    make_collapse_map(old->nsamples, map, collapse);    
    
    // Create a new BranchTab
    BranchTab *new = BranchTab_new(newsamp);

    for(unsigned i=1; i < old->dim; ++i) {
        tipId_t tid = remap_bits(old->nsamples, map, i);
        BranchTab_add(new, tid, old->tab[i]);
    }
    return new;
}

/// Remove populations
BranchTab *BranchTab_rmPops(BranchTab *old, tipId_t remove) {
    int newsamp = old->nsamples - num1bits(remove);

    // Make map, an array whose i'th entry is an unsigned integer
    // with one bit on and the rest off. The on bit indicates the
    // position in the new id of the i'th bit in the old id.
    tipId_t map[old->nsamples];
    make_rm_map(old->nsamples, map, remove);
    
    // Create a new BranchTab
    BranchTab *new = BranchTab_new(newsamp);
    CHECKMEM(new);
    for(unsigned i=1; i < old->dim; ++i) {
        tipId_t tid = remap_bits(old->nsamples, map, i);
        if(id)
            BranchTab_add(new, tid, old->tab[i]);
    }
    return new;
}

/// Return sum of values in BranchTab.
long double BranchTab_sum(const BranchTab *self) {
    long double s=0;

    for(unsigned i = 1; i < self->dim; ++i)
        s += self->tab[i];

    return s;
}

/// Return negative of sum of p*ln(p)
long double BranchTab_entropy(const BranchTab *self) {
    assert(Dbl_near(1.0, BranchTab_sum(self)));
    long double entropy=0.0;

    for(unsigned i = 1; i < self->dim; ++i) {
        long double x = self->tab[i];
        entropy -= x * logl(x);
    }

    return entropy;
}

/// Divide all values by their sum. Return 0 on success, or 1 on
/// failure.
int BranchTab_normalize(BranchTab *self) {
    long double s = BranchTab_sum(self);

    if(s==0)
        return 1;

    // divide by sum
    for(unsigned i = 1; i < self->dim; ++i) {
        self->tab[i] /= s;
    }

    return 0;
}

/// Calculate KL divergence from two BranchTab objects, which
/// should be normalized before entering this function. Use
/// BranchTab_normalize to normalize. Function returns HUGE_VAL if there
/// are observed values without corresponding values in expt.
long double BranchTab_KLdiverg(const BranchTab *obs, const BranchTab *expt) {
    assert(LDbl_near(1.0L, BranchTab_sum(obs)));
    assert(LDbl_near(1.0L, BranchTab_sum(expt)));
    assert(obs->dim == expt->dim);

    long double kl=0.0;
    long double p;  // observed frequency
    long double q;  // frequency under model
    for(unsigned i=1; i < obs->dim; ++i) {
        p = obs->tab[i];
        q = expt->tab[i];

        // Use p and q to add a term to kl.
        if(p == 0.0) {
            // Do nothing: p*log(p/q) -> 0 as p->0, regardless of
            // q. This is because p*log(p/q) is the log of
            // (p/q)**p, which equals 1 if p=0, no matter the value
            // of q.
        }else{
            if(q==0.0)
                return HUGE_VAL;
            kl += p*logl(p/q);
        }
    }
    return kl;
}

/// Negative log likelihood. Multinomial model.
/// lnL is sum across site patterns of x*log(p), where x is an
/// observed site pattern count and p its probability.
long double BranchTab_negLnL(const BranchTab *obs, const BranchTab *expt) {
    assert(LDbl_near(1.0L, BranchTab_sum(expt)));
    assert(obs->dim == expt->dim);

    long double lnL=0.0;
    long double x;  // observed count
    long double p;  // probability under model
    for(unsigned i=0; i < obs->dim; ++i) {
        x = obs->tab[i];
        p = expt->tab[i];
        if(p == 0.0) {
            if(x != 0.0)
                return HUGE_VAL;  // blows up
        }else
            lnL += x*logl(p);
    }
    return -lnL;
}

/// Make map, an array whose i'th entry is an unsigned integer with one
/// bit on and the rest off. The on bit indicates the position in the
/// new id of the i'th bit in the old id.
///
/// All bits equal to 1 in collapse are mapped to the minimum bit in
/// collapse. All bits with positions below that of this mininum bit
/// are mapped to their original position. Bits that are above the
/// minumum bit and are off in collapse are shifted right by "shift"
/// places, where shift is one less than the number of on bits in
/// collapse that are to the right of the bit in question.
///
/// On entry, n equals the number of samples, map is an array of n
/// tipId_t values, and collapse is an integer whose "on" bits
/// represent the set of samples that will be collapsed into a single
/// sample. 
void make_collapse_map(size_t n, tipId_t map[n], tipId_t collapse) {
    unsigned maxbit = 8*sizeof(tipId_t) - nlz(collapse);
    if(maxbit > n) {
        fprintf(stderr,"%s:%s:%d: maxbit (=%u) cannot exceed n (=%zu)\n",
                __FILE__,__func__,__LINE__, maxbit, n);
        exit(EXIT_FAILURE);
    }
    int i, shift=0;
    tipId_t min = 0, bit = 1u;
    for(i=0; i < n; ++i, bit <<= 1) {
        if( collapse & bit ) {
            if(min == 0) {
                min = bit;
            }else
                ++shift;
            map[i] = min;
        }else
            map[i] = bit >> shift;
    }
}

/// Make map, an array whose i'th entry is an unsigned integer with one
/// bit on and the rest off. The on bit indicates the position in the
/// new id of the i'th bit in the old id.
///
/// All samples corresponding to an "on" bit in "remove" are deleted.
/// Other samples (those not removed) have their bits right-shifted
/// by a number of positions equal to the number of lower-order
/// populations that have been removed.
///
/// On entry, n equals the number of samples, map is an array of n
/// tipId_t values, and remove is an integer whose "on" bits
/// represent the set of samples that will be removed.
void make_rm_map(size_t n, tipId_t map[n], tipId_t remove) {
    unsigned maxbit = 8*sizeof(tipId_t) - nlz(remove);
    if(maxbit > n) {
        fprintf(stderr,"%s:%s:%d: maxbit (=%u) cannot exceed n (=%zu)\n",
                __FILE__,__func__,__LINE__, maxbit, n);
        exit(EXIT_FAILURE);
    }
    int i, shift=0;
    tipId_t bit = 1u;
    for(i=0; i < n; ++i, bit <<= 1) {
        if( remove & bit ) {
            ++shift;
            map[i] = 0;
        }else
            map[i] = bit >> shift;
    }
}

// Remap the bits in a tipId_t variable. Array "map" specifies how
// the bits should be rearranged. Function returns the remapped
// value of "old".
tipId_t remap_bits(size_t n, tipId_t map[n], tipId_t old) {
    tipId_t new = 0u, bit=1u;
    for(int i=0; i<n; ++i, bit <<= 1)
        if( old & bit )
            new |= map[i];
    return new;
}

#ifdef TEST

#include "binary.h"
#include <string.h>
#include <assert.h>
#include <unistd.h>

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

    unsigned nsamples = 5;
    unsigned maxtid = low_bits_on(nsamples) - 1;
    BranchTab *bt = BranchTab_new(nsamples);
    CHECKMEM(bt);
    assert(low_bits_on(nsamples) == BranchTab_size(bt));

    tipId_t key[maxtid];
    long double  val[maxtid];

    int i;
    for(i=1; i <= maxtid ; ++i) {
        key[i-1] = i;
        val[i-1] = (long double) i;
        BranchTab_add(bt, key[i-1], val[i-1]);
    }
    assert(BranchTab_hasSingletons(bt));

    for(i=0; i < maxtid; ++i) {
        assert(val[i] == BranchTab_get(bt, key[i]));
    }

    long double x = BranchTab_sum(bt);
    if(verbose)
        printf("BranchTab_sum: %Lg\n", x);
    assert(x == maxtid*((maxtid+1)/2.0L));
    assert(0 == BranchTab_normalize(bt));
    assert(LDbl_near(1.0, BranchTab_sum(bt)));
    x = BranchTab_entropy(bt);
    if(verbose)
        printf("BranchTab_entropy: %Lg\n", x);

    if(verbose)
        BranchTab_print(bt, stdout);

    BranchTab *bt2 = BranchTab_dup(bt);
    assert(BranchTab_equals(bt, bt2));

    if(verbose)
        BranchTab_print(bt, stdout);

    // Test  _plusEquals, _minusEquals, and _dup
    BranchTab_minusEquals(bt, bt2);
    for(i=0; i<maxtid; ++i)
        assert(0.0 == BranchTab_get(bt, key[i]));

    BranchTab_free(bt2);
    bt2 = BranchTab_collapse(bt, 03);
    if(verbose) {
        printf("size before collapse: %u; after: %u\n",
               BranchTab_size(bt), BranchTab_size(bt2));
        BranchTab_print(bt2, stdout);
    }
    assert(BranchTab_size(bt) > BranchTab_size(bt2));
    assert(LDbl_near(BranchTab_sum(bt), BranchTab_sum(bt2)));
    unitTstResult("BranchTab_collapse", "OK");
    
    // test make_collapse_map and remap_bits
    size_t n = 8*sizeof(tipId_t);
    tipId_t map[n];
    tipId_t id1, id2, id3;
    id1 = 04444;
    id2 = 033333;
    make_collapse_map(n, map, id1);
    id3 = remap_bits(n, map, id2);
    if(verbose) {
        printf("After make_collapse_map and remap_bits...\n");
        puts("id1: ");
        printBits(sizeof(id1), &id1, stdout);
        puts("id2: ");
        printBits(sizeof(id2), &id2, stdout);
        puts("id3: ");
        printBits(sizeof(id3), &id3, stdout);
    }
    assert(04 == remap_bits(n, map, id1));
    assert(04 == remap_bits(n, map, 040));
    assert(04 == remap_bits(n, map, 0400));
    assert(04 == remap_bits(n, map, 04000));
    assert(0100 == remap_bits(n, map, 0200));
    assert(0400 == remap_bits(n, map, 02000));
    assert(02000 == remap_bits(n, map, 020000));
    assert(03773 == remap_bits(n, map, id2));
    unitTstResult("make_collapse_map", "OK");
    unitTstResult("remap_bits", "OK");

    // test make_rm_map
    memset(map, 0, sizeof(map));
    id1 = 044;
    id2 = 077;
    make_rm_map(n, map, id1);
    id3 = remap_bits(n, map, id2);
    if(verbose) {
        printf("After make_rm_map and remap_bits...\n");
        puts("id1: ");
        printBits(sizeof(id1), &id1, stdout);
        puts("id2: ");
        printBits(sizeof(id2), &id2, stdout);
        puts("id3: ");
        printBits(sizeof(id3), &id3, stdout);
    }
    assert(017 == id3);
    assert(03 == remap_bits(n, map, 07));
    assert(04 == remap_bits(n, map, 010));
    unitTstResult("make_rm_map", "OK");

    BranchTab_free(bt2);
    bt2 = BranchTab_rmPops(bt, 06);
    if(verbose) {
        printf("size before rmPops: %u; after: %u\n",
               BranchTab_size(bt), BranchTab_size(bt2));
        BranchTab_print(bt2, stdout);
    }

    assert(BranchTab_size(bt) > BranchTab_size(bt2));
    assert(BranchTab_sum(bt) >= BranchTab_sum(bt2));
    unitTstResult("BranchTab_rmPops", "OK");

    BranchTab_free(bt);
    BranchTab_free(bt2);
    unitTstResult("BranchTab", "OK");

    return 0;
}

#endif
