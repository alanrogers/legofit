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
#include "tokenizer.h"
#include "lblndx.h"
#include "parstore.h"
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/// Dimension of hash table. Must be a power of 2
#define BT_DIM 128u

/// Make sure BT_DIM is a power of 2
#if (BT_DIM==0u || (BT_DIM & (BT_DIM-1u)))
# error BT_DIM must be a power of 2
#endif

/// A single element of hash table
typedef struct BTLink {
    struct BTLink  *next;
    tipId_t         key;
    double          value;
} BTLink;

/// Hash table for branch lengths.
struct BranchTab {
#ifndef NDEBUG
    int frozen;   // nonzero => no further changes allowed
#endif
    BTLink         *tab[BT_DIM];
};

BTLink     *BTLink_new(tipId_t key, double value);
BTLink     *BTLink_add(BTLink * self, tipId_t key, double value);
double      BTLink_get(BTLink * self, tipId_t key);
int         BTLink_hasSingletons(BTLink * self);
void        BTLink_free(BTLink * self);
void        BTLink_printShallow(const BTLink * self, FILE *fp);
void        BTLink_print(const BTLink * self, FILE *fp);
BTLink     *BTLink_dup(const BTLink *self);
int         BTLink_equals(const BTLink *lhs, const BTLink *rhs);
static void make_map(size_t n, tipId_t map[n], tipId_t collapse);
static void make_rm_map(size_t n, tipId_t map[n], tipId_t remove);
static tipId_t remap_bits(size_t n, tipId_t map[n], tipId_t old);

#if TIPID_SIZE==32
static inline uint32_t tipIdHash(uint32_t key);
#elif TIPID_SIZE==64
static inline uint32_t tipIdHash(uint64_t key);
#endif

#if TIPID_SIZE==32
static inline uint32_t tipIdHash( uint32_t key) {
    return uint32Hash(key);
}
#elif TIPID_SIZE==64
static inline uint64_t tipIdHash( uint64_t key) {
    return uint64Hash(key);
}
#else
#error "Can't compile tipIdHash function. See branchtab.c"
#endif

/// Constructor of class BTLink
BTLink         *BTLink_new(tipId_t key, double value) {
    BTLink         *new = malloc(sizeof(*new));
    CHECKMEM(new);

    new->next = NULL;
    new->key = key;
    new->value = value;
    return new;
}

/// Duplicate linked list.
BTLink     *BTLink_dup(const BTLink *old) {
    if(old == NULL)
        return NULL;

    BTLink *new = malloc(sizeof(*new));
    CHECKMEM(new);
    new->key = old->key;
    new->value = old->value;
    new->next = BTLink_dup(old->next);
    return new;
}

/// Return 1 if the two linked lists are equal, or 0 otherwise.
int  BTLink_equals(const BTLink *lhs, const BTLink *rhs) {
    if(lhs==NULL && rhs==NULL)
        return 1;
    if(lhs==NULL || rhs==NULL)
        return 0;
    if(lhs->key!=rhs->key
       || lhs->value!=rhs->value)
        return 0;
    return BTLink_equals(lhs->next, rhs->next);
}

/// Destructor for BTLink
void BTLink_free(BTLink * self) {
    if(self == NULL)
        return;
    BTLink_free(self->next);
    free(self);
}

/// Add a value to a BTLink object. On return, self->value
/// equals the old value plus the new one.
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

/// Return 1 if linked list contains keys for singletone site patterns.
int BTLink_hasSingletons(BTLink * self) {
    if(self == NULL)
        return 0;
    if(isPow2(self->key))
        return 1;
    return BTLink_hasSingletons(self->next);
}

void BTLink_printShallow(const BTLink * self, FILE *fp) {
    if(self == NULL)
        return;
    fprintf(fp, " [%lo, %lf]", (unsigned long) self->key, self->value);
}

void BTLink_print(const BTLink * self, FILE *fp) {
    if(self == NULL)
        return;
    BTLink_printShallow(self, fp);
    BTLink_print(self->next, fp);
}

BranchTab    *BranchTab_new(void) {
    BranchTab    *new = malloc(sizeof(*new));
    CHECKMEM(new);
    memset(new, 0, sizeof(*new));
    return new;
}

BranchTab    *BranchTab_dup(const BranchTab *old) {
    BranchTab *new = BranchTab_new();

#ifndef NDEBUG
    new->frozen = old->frozen;
#endif

    int i;
    for(i=0; i < BT_DIM; ++i)
        new->tab[i] = BTLink_dup(old->tab[i]);

    return new;
}

/// Return 1 if two BranchTab objects are equal; 0 otherwise.
int BranchTab_equals(const BranchTab *lhs, const BranchTab *rhs) {
    int i;

    for(i=0; i < BT_DIM; ++i) {
        if(!BTLink_equals(lhs->tab[i], rhs->tab[i]))
            return 0;
    }
#ifndef NDEBUG
    if(lhs->frozen != rhs->frozen)
        return 0;
#endif
    return 1;
}

/// Destructor for BranchTab.
void BranchTab_free(BranchTab * self) {
    int i;
    for(i=0; i < BT_DIM; ++i)
        BTLink_free(self->tab[i]);
    free(self);
}

/// Return 1 if BranchTab includes singleton site patterns
int BranchTab_hasSingletons(BranchTab * self) {
    int i;
    for(i=0; i < BT_DIM; ++i) {
        if(BTLink_hasSingletons(self->tab[i]))
            return 1;
    }
    return 0;
}

/// Return value corresponding to key, or nan if no value is found.
double BranchTab_get(BranchTab * self, tipId_t key) {
    unsigned h = tipIdHash(key) & (BT_DIM-1u);
    assert(h < BT_DIM);
    assert(self);
    return BTLink_get(self->tab[h], key);
}

/// Add a value to table. If key already exists, new value is added to
/// old one.
void BranchTab_add(BranchTab * self, tipId_t key, double value) {
    assert(!self->frozen);
    if(key==0) {
        fprintf(stderr,"%s:%d: tipId_t mustn't equal 0\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    unsigned h = tipIdHash(key) & (BT_DIM-1u);
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

/// Divide all values by denom. Return 0 on success, or 1 on failure.
int BranchTab_divideBy(BranchTab *self, double denom) {
#ifndef NDEBUG
    assert(!self->frozen);
    self->frozen = 1;  // you can only call this function once
#endif

    // divide by denom
    unsigned i;
    for(i = 0; i < BT_DIM; ++i) {
        BTLink *el;
        for(el = self->tab[i]; el; el = el->next)
            el->value /= denom;
    }

    return 0;
}

/// Print a BranchTab to standard output.
void BranchTab_print(const BranchTab *self, FILE *fp) {
    unsigned i;
    for(i=0; i < BT_DIM; ++i) {
        fprintf(fp, "%2u:", i);
        BTLink_print(self->tab[i], fp);
        putc('\n', fp);
    }
}

/// Add each entry in table rhs to table lhs
void BranchTab_plusEquals(BranchTab *lhs, BranchTab *rhs) {
    assert(!lhs->frozen && !rhs->frozen);
    int i;
    for(i=0; i<BT_DIM; ++i) {
        BTLink *link;
        for(link=rhs->tab[i]; link!=NULL; link=link->next)
            BranchTab_add(lhs, link->key, link->value);
    }
}

/// Subtract each entry in table rhs to table lhs
void BranchTab_minusEquals(BranchTab *lhs, BranchTab *rhs) {
    assert(!lhs->frozen && !rhs->frozen);
    int i;
    for(i=0; i<BT_DIM; ++i) {
        BTLink *link;
        for(link=rhs->tab[i]; link!=NULL; link=link->next) {
            double x = BranchTab_get(lhs, link->key);
            if(isnan(x)) {
                fprintf(stderr,"%s:%s:%d: incompatible arguments\n",
                        __FILE__,__func__,__LINE__);
                exit(EXIT_FAILURE);
            }
            BranchTab_add(lhs, link->key, -link->value);
        }
    }
}

/// Fill arrays key, value, and square with values in BranchTab.
/// On return, key[i] is the id of the i'th site pattern, value[i] is
/// the total branch length associated with that site pattern, and
/// sumsqr[i] is the corresponding sum of squared branch lengths.
void BranchTab_toArrays(BranchTab *self, unsigned n, tipId_t key[n],
                        double value[n], double sumsqr[n]) {
    int i, j=0;
    for(i=0; i<BT_DIM; ++i) {
        BTLink *link;
        for(link=self->tab[i]; link!=NULL; link=link->next) {
            if(j >= n)
                eprintf("%s:%s:%d: buffer overflow\n",
                        __FILE__,__func__,__LINE__);
            key[j] = link->key;
            value[j] = link->value;
            ++j;
        }
    }
}

/// Construct a BranchTab by parsing an input file.
/// Recognizes comments, which extend from '#' to end-of-line.
BranchTab *BranchTab_parse(const char *fname, const LblNdx *lblndx) {
    FILE *fp = efopen(fname, "r");

    BranchTab *self = BranchTab_new();
    CHECKMEM(self);

    int         ntokens;
    char        buff[500];
    Tokenizer  *tkz = Tokenizer_new(50);

    while(1) {
        if(fgets(buff, sizeof(buff), fp) == NULL)
            break;

        if(!strchr(buff, '\n') && !feof(fp))
            eprintf("s:%s:%d: buffer overflow. buff size: %zu\n",
                    __FILE__, __func__, __LINE__, sizeof(buff));

        // strip trailing comments
        char *comment = strchr(buff, '#');
        if(comment)
            *comment = '\0';

        Tokenizer_split(tkz, buff, " \t"); // tokenize
        ntokens = Tokenizer_strip(tkz, " \t\n");
        if(ntokens == 0)
            continue;

        char *tok = Tokenizer_token(tkz, 0);
        tipId_t key=LblNdx_getTipId(lblndx, tok);
        if(key==0) {
            fprintf(stderr,"%s:%s:%d: can't find id for label %s\n",
                    __FILE__,__func__,__LINE__, tok);
            exit(EXIT_FAILURE);
        }

        tok = Tokenizer_token(tkz, 1);
        errno = 0;
        double prob = strtod(tok, NULL);
        if(errno) {
            fprintf(stderr,"%s:%s:%d: Can't parse 2nd field as float.\n",
                    __FILE__,__func__,__LINE__);
            fprintf(stderr," input:");
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }

        BranchTab_add(self, key, prob);
    }
    fclose(fp);
    Tokenizer_free(tkz);
    return self;
}

/// Map two or more populations into a single population.
BranchTab *BranchTab_collapse(BranchTab *old, tipId_t collapse) {
    int n = 8 * sizeof(tipId_t); // number of bits

    // Make map, an array whose i'th entry is an unsigned integer
    // with one bit on and the rest off. The on bit indicates the
    // position in the new id of the i'th bit in the old id.
    tipId_t map[n];
    make_map(n, map, collapse);

    // Create a new BranchTab
    BranchTab *new = BranchTab_new();
    for(int i=0; i < BT_DIM; ++i) {
        BTLink *el;
        for(el = old->tab[i]; el; el = el->next)
            BranchTab_add(new, remap_bits(n, map, el->key), el->value);
    }
    return new;
}

/// Remove populations
BranchTab *BranchTab_rmPops(BranchTab *old, tipId_t remove) {
    int n = 8 * sizeof(tipId_t); // number of bits

    // Make map, an array whose i'th entry is an unsigned integer
    // with one bit on and the rest off. The on bit indicates the
    // position in the new id of the i'th bit in the old id.
    tipId_t map[n];
    make_rm_map(n, map, remove);

    // Create a new BranchTab
    BranchTab *new = BranchTab_new();
    for(int i=0; i < BT_DIM; ++i) {
        BTLink *el;
        for(el = old->tab[i]; el; el = el->next) {
            tipId_t id = remap_bits(n, map, el->key);
            if(id)
                BranchTab_add(new, id, el->value);
        }
    }
    return new;
}

// Make map, an array whose i'th entry is an unsigned integer with one
// bit on and the rest off. The on bit indicates the position in the
// new id of the i'th bit in the old id.
//
// All bits equal to 1 in collapse are mapped to the minimum bit in
// collapse. All bits with positions below that of this mininum bit
// are mapped to their original position. Aits that are above the
// minumum bit and are off in collapse are shifted right by "shift"
// places, where shift is one less than the number of on bits in
// collapse that are to the right of the bit in question.
static void make_map(size_t n, tipId_t map[n], tipId_t collapse) {
    int i, shift=0;
    tipId_t min = n, bit = 1u;
    for(i=0; i < n; ++i, bit <<= 1) {
        if( collapse & bit ) {
            if(min == n) {
                min = bit;
            }else
                ++shift;
            map[i] = min;
        }else
            map[i] = bit >> shift;
    }
}

static void make_rm_map(size_t n, tipId_t map[n], tipId_t remove) {
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
static tipId_t remap_bits(size_t n, tipId_t map[n], tipId_t old) {
    assert(n == 8*sizeof(tipId_t));
    tipId_t new = 0u, bit=1u;
    for(int i=0; i<n; ++i, bit <<= 1)
        if( old & bit )
            new |= map[i];
    return new;
}

/// Return sum of values in BranchTab.
double BranchTab_sum(const BranchTab *self) {
    unsigned i;
    double s=0.0;

    for(i = 0; i < BT_DIM; ++i) {
        BTLink *el;
        for(el = self->tab[i]; el; el = el->next)
            s += el->value;
    }

    return s;
}

/// Return negative of sum of p*ln(p)
double BranchTab_entropy(const BranchTab *self) {
    assert(Dbl_near(1.0, BranchTab_sum(self)));
    unsigned i;
    double entropy=0.0;

    for(i = 0; i < BT_DIM; ++i) {
        BTLink *el;
        for(el = self->tab[i]; el; el = el->next)
            entropy -= el->value * log(el->value);
    }

    return entropy;
}

/// Divide all values by their sum. Return 0
/// on success, or 1 on failure.
int BranchTab_normalize(BranchTab *self) {
    unsigned i;
    double s = BranchTab_sum(self);

    if(s==0)
        return 1;

    // divide by sum
    for(i = 0; i < BT_DIM; ++i) {
        BTLink *el;
        for(el = self->tab[i]; el; el = el->next)
            el->value /= s;
    }

    return 0;
}

/// Calculate KL divergence from two BranchTab objects, which
/// should be normalized before entering this function. Use
/// BranchTab_normalize to normalize. Function returns HUGE_VAL if there
/// are observed values without corresponding values in expt.
double BranchTab_KLdiverg(const BranchTab *obs, const BranchTab *expt) {
    assert(Dbl_near(1.0, BranchTab_sum(obs)));
    assert(Dbl_near(1.0, BranchTab_sum(expt)));

    int i;
    double kl=0.0;
    double p;  // observed frequency
    double q;  // frequency under model
    for(i=0; i < BT_DIM; ++i) {
        BTLink *o, *e;
        o = obs->tab[i];
        e = expt->tab[i];
        while(o && e) {
            if(o->key < e->key) {
                p = o->value;
                q = 0.0;
                o = o->next;
            }else if(o->key > e->key) {
                p = 0.0;
                q = e->value;
                e = e->next;
            }else {
                assert(o->key == e->key);
                p = o->value;
                q = e->value;
                e = e->next;
                o = o->next;
            }
            if(p == 0.0) {
                // Do nothing: p*log(p/q) -> 0 as p->0, regardless of
                // q. This is because p*log(p/q) is the log of
                // (p/q)**p, which equals 1 if p=0, no matter the value
                // of q.
            }else{
                if(q==0.0)
                    return HUGE_VAL;
                kl += p*log(p/q);
            }
        }
        while(o) { // e->value is q=0: fail unless p=0
            p = o->value;
            if(p != 0.0)
                return HUGE_VAL;
            o = o->next;
        }
        // Any remaining cases have p=0, so contribution to kl is 0.
    }
    return kl;
}

/// Negative log likelihood. Multinomial model.
/// lnL is sum across site patterns of x*log(p), where x is an
/// observed site pattern count and p its probability.
double BranchTab_negLnL(const BranchTab *obs, const BranchTab *expt) {
    assert(Dbl_near(1.0, BranchTab_sum(expt)));

    int i;
    double lnL=0.0;
    double x;  // observed count
    double p;  // probability under model
    for(i=0; i < BT_DIM; ++i) {
        BTLink *o, *e;  // observed and expected
        o = obs->tab[i];
        e = expt->tab[i];
        while(o && e) {
            if(o->key < e->key) {
                x = o->value;
                p = 0.0;
                o = o->next;
            }else if(o->key > e->key) {
                x = 0.0;
                p = e->value;
                e = e->next;
            }else {
                assert(o->key == e->key);
                x = o->value;
                p = e->value;
                e = e->next;
                o = o->next;
            }
            if(p == 0.0) {
                if(x != 0.0)
                    return HUGE_VAL;  // blows up
            }else
                lnL += x*log(p);
        }
        while(o) { // e->value is p=0
            x = o->value;
            if(x != 0.0)
                return HUGE_VAL; // blows up
            o = o->next;
        }
        // Any remaining cases have x=0, so contribution to lnL is 0.
    }
    return -lnL;
}


#ifdef TEST

#include "binary.h"
#include <string.h>
#include <assert.h>
#include <unistd.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

//      a-------|
//              |ab--|
//      b--|bb--|    |
//         |         |abc--
//         |c--------|
//
//  t = 0  1    3    5.5     inf
const char *tstInput =
    " # this is a comment\n"
    "time fixed  T0=0\n"
    "time free   Tc=1\n"
    "time free   Tab=3\n"
    "time free   Tabc=5.5\n"
    "twoN free   twoNa=100\n"
    "twoN fixed  twoNb=123\n"
    "twoN free   twoNc=213.4\n"
    "twoN fixed  twoNbb=32.1\n"
    "twoN free   twoNab=222\n"
    "twoN fixed  twoNabc=1.2e2\n"
    "mixFrac free Mc=0.02\n"
    "segment a   t=T0     twoN=twoNa    samples=1\n"
    "segment b   t=T0     twoN=twoNb    samples=1\n"
    "segment c   t=Tc     twoN=twoNc    samples=1\n"
    "segment bb  t=Tc     twoN=twoNbb\n"
    "segment ab  t=Tab    twoN=twoNab\n"
    "segment abc t=Tabc   twoN=twoNabc\n"
    "mix    b  from bb + Mc * c\n"
    "derive a  from ab\n"
    "derive bb from ab\n"
    "derive ab from abc\n"
    "derive c  from abc\n";
const char *tstPatProbInput =
    "#SitePat   obs\n"
    "a:b        2.0\n"
    "a:c        1.0\n"
    "b:c        1.0\n";

int main(int argc, char **argv) {
    int verbose=0;
    if(argc > 1) {
        if(argc!=2 || 0!=strcmp(argv[1], "-v")) {
            fprintf(stderr,"usage: xbranchtab [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    const char *tstFname = "mktree-tmp.lgo";
    FILE       *fp = fopen(tstFname, "w");
    fputs(tstInput, fp);
    fclose(fp);

    const char *tstPatProbFname = "patprob-tmp.txt";
    fp = fopen(tstPatProbFname, "w");
    fputs(tstPatProbInput, fp);
    fclose(fp);

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

    double x = BranchTab_sum(bt);
    if(verbose)
        printf("BranchTab_sum: %lg\n", x);
    assert(0 == BranchTab_normalize(bt));
    assert(Dbl_near(1.0, BranchTab_sum(bt)));
    x = BranchTab_entropy(bt);
    if(verbose)
        printf("BranchTab_entropy: %lg\n", x);

    if(verbose)
        BranchTab_print(bt, stdout);


    BranchTab *bt2 = BranchTab_dup(bt);
    assert(BranchTab_equals(bt, bt2));

    if(verbose)
        BranchTab_print(bt, stdout);

    // Test  _plusEquals, _minusEquals, and _dup
    BranchTab_minusEquals(bt, bt2);
    for(i=0; i<25; ++i)
        assert(0.0 == BranchTab_get(bt, key[i]));

    BranchTab_free(bt);

    Bounds   bnd = {
                    .lo_twoN = 0.0,
                    .hi_twoN = 1e7,
                    .lo_t = 0.0,
                    .hi_t = HUGE_VAL
    };
    GPTree *g = GPTree_new(tstFname, bnd);
    LblNdx lblndx = GPTree_getLblNdx(g);

    bt = BranchTab_parse(tstPatProbFname, &lblndx);

    // test make_map and remap_bits
    size_t n = 8*sizeof(tipId_t);
    tipId_t map[n];
    tipId_t id1, id2, id3;
    id1 = 04444;
    id2 = 033333;
    make_map(n, map, id1);
    id3 = remap_bits(n, map, id2);
    if(verbose) {
        printf("After make_map and remap_bits...\n");
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
    unitTstResult("make_map", "OK");
    unitTstResult("remap_bits", "OK");

    BranchTab_free(bt2);
    bt2 = BranchTab_collapse(bt, 03);
    if(verbose) {
        printf("size before collapse: %u; after: %u\n",
               BranchTab_size(bt), BranchTab_size(bt2));
        BranchTab_print(bt2, stdout);
    }
    assert(BranchTab_size(bt) > BranchTab_size(bt2));
    assert(Dbl_near(BranchTab_sum(bt), BranchTab_sum(bt2)));
    unitTstResult("BranchTab_collapse", "OK");
    
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
    GPTree_free(g);
    unitTstResult("BranchTab", "untested");
    unlink(tstFname);
    unlink(tstPatProbFname);
}
#endif
