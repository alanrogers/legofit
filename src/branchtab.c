/**
 * @file branchtab.c
 * @author Alan R. Rogers
 * @brief Hash table associating key (an unsigned int whose bits encode
 * the descendants of a node in a gene tree) and value (a double
 * representing the length of the ascending branch)
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
#if COST==CHISQR_COST
    double          sqr;   // squares
#endif
} BTLink;

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
void        BTLink_printShallow(const BTLink * self);
void        BTLink_print(const BTLink * self);
BTLink     *BTLink_dup(const BTLink *self);
int         BTLink_equals(const BTLink *lhs, const BTLink *rhs);

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
#if COST==CHISQR_COST
    new->sqr = value*value;
#endif
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
#if COST==CHISQR_COST
    new->sqr = old->sqr;
#endif
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
#if COST==CHISQR_COST
    if(lhs->sqr != rhs->sqr)
        return 0;
#endif
    return BTLink_equals(lhs->next, rhs->next);
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
#if COST==CHISQR_COST
    self->sqr += value*value;
#endif
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

void BTLink_printShallow(const BTLink * self) {
    if(self == NULL)
        return;
#if COST==CHISQR_COST
    printf(" [%lu, %lf, %lf]",
           (unsigned long) self->key, self->value, self->sqr);
#else
    printf(" [%lu, %lf]", (unsigned long) self->key, self->value);
#endif
}

void BTLink_print(const BTLink * self) {
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
    unsigned h = tipIdHash(key);
    assert(h < BT_DIM);
    assert(self);
    return BTLink_get(self->tab[h], key);
}

/// Add a value to table. If key already exists, new value is added to
/// old one.  
void BranchTab_add(BranchTab * self, tipId_t key, double value) {
    assert(!self->frozen);
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
        for(el = self->tab[i]; el; el = el->next) {
            el->value /= denom;
#if COST==CHISQR_COST
            el->sqr /= denom;
#endif
        }
    }

    return 0;
}

void BranchTab_print(const BranchTab *self) {
    unsigned i;
    for(i=0; i < BT_DIM; ++i) {
        printf("%2u:", i);
        BTLink_print(self->tab[i]);
        putchar('\n');
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

void BranchTab_toArrays(BranchTab *self, unsigned n, tipId_t key[n],
						double value[n], double sqr[n]) {
    int i, j=0;
    for(i=0; i<BT_DIM; ++i) {
        BTLink *link;
        for(link=self->tab[i]; link!=NULL; link=link->next) {
            if(j >= n)
                eprintf("%s:%s:%d: buffer overflow\n",
						__FILE__,__func__,__LINE__);
            key[j] = link->key;
            value[j] = link->value;
#if COST==CHISQR_COST
            sqr[j] = link->sqr;
#endif
            ++j;
        }
    }
}

/// Construct by parsing an input file. Returns a pointer
/// to a normalized BranchTab object.
BranchTab *BranchTab_parse(const char *fname, const LblNdx *lblndx) {
    FILE *fp = fopen(fname, "r");
    if(fp == NULL)
        eprintf("%s:%d: Can't read file \"%s\".\n",
                __FILE__, __LINE__, fname);

    BranchTab *self = BranchTab_new();
    CHECKMEM(self);

    int         i, ntokens;
    char        buff[500];
    char        lblbuff[100];
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

        char *tok = Tokenizer_token(tkz, 1);
        errno = 0;
        double prob = strtod(tok, NULL);
        if(errno) {
            fprintf(stderr,"%s:%s:%d: Can't parse 2nd field as float.\n",
                    __FILE__,__func__,__LINE__);
            fprintf(stderr," input:");
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }

        snprintf(lblbuff, sizeof(lblbuff), "%s", Tokenizer_token(tkz, 0));
        Tokenizer_split(tkz, lblbuff, ":");
        ntokens = Tokenizer_strip(tkz, " \t\n");
        assert(ntokens > 0);
        tipId_t key=0, id;

        // Get tipId corresponding to label
        for(i=0; i < ntokens; ++i) {
            tok = Tokenizer_token(tkz, i);
            id = LblNdx_getTipId(lblndx, tok);
            if(id == 0)
                eprintf("%s:%s:%d: unrecognized label, %s, in input.\n",
                        __FILE__,__func__,__LINE__, tok);
            key |= id;
        }
        BranchTab_add(self, key, prob);
    }
    fclose(fp);
    Tokenizer_free(tkz);
    return self;
}

#if COST==CHISQR_COST
/// Chi-squared difference between observed and expected.  The
/// Chi-squared statistic is (obs - expected)^2/variance. To estimate
/// the variance, note that the expected number of mutations is
/// u*nnuc*Poisson(B), where B is a random variable, the branch length per
/// nucleotide site contributing to a given site pattern. This is u*nnuc
/// times a mixture of Poisson distributions. It's mean is u*nnuc*Mean(B),
/// and its variance is u*nnuc*Mean(B) + nnuc*u*u*Var(B).
double        BranchTab_cost(const BranchTab *obs, const BranchTab *expt,
                             double u, long nnuc, double n) {
    assert(expt->frozen);
    int i;
    double U = u*nnuc;
    double cost=0.0, diff;
    double v;     // nnuc*u*u*Var(B), the variance of branch length
    double obval;  // observed mutations
    double exval;  // mutations expected under model
    for(i=0; i < BT_DIM; ++i) {
        BTLink *o = obs->tab[i];
        BTLink *e = expt->tab[i];
        while(o && e) {
            if(o->key < e->key) {
                // e link missing, so exval=0.
                obval = o->value;
                exval = 0.0;
                o = o->next;
                v = 0.0;
            }else if(o->key > e->key) {
                // o link missing, so obval=0.
                obval = 0.0;
                exval = e->value * U;
                v = e->sqr - e->value * e->value;
                assert(v>=0.0);
                v *= u*U*n/(n-1.0);
                e = e->next;
            }else {
                assert(o->key == e->key);
                obval = o->value;
                exval = e->value * U;
                v = e->sqr - e->value * e->value;
                assert(v>=0.0);
                v *= u*U*n/(n-1.0);
                e = e->next;
                o = o->next;
            }
            diff = obval-exval;
#if 0
            printf("o=%lg e=%lg v=%lg chisq=%lg\n",
                   obval, exval, v, diff*diff/(exval+v));
#endif
            cost += diff*diff/(exval+v);
        }
        if(o) { // e link missing, so exval=0
            cost = HUGE_VAL;
            break;
        }
        obval=0.0;
        while(e) { // o link missing, so obval=0
            exval = e->value * U;
            v = e->sqr - e->value * e->value;
            assert(v>=0.0);
            v *= u*U*n/(n-1.0);
            diff = obval-exval;
#if 0
            printf("o=%lg e=%lg v=%lg chisq=%lg\n",
                   obval, exval, v, diff*diff/(exval+v));
#endif
            assert(v>=0.0);
            cost += diff*diff/(exval+v);
            e = e->next;
        }
    }
    assert(cost >= 0.0);
    return cost;
}
#endif

#if COST==POISSON_COST
/// Use Poisson model to calculate negative log likelihood. 
double BranchTab_poissonCost(const BranchTab *obs, const BranchTab *expt,
                             double u, long nnuc, double n) {
    assert(expt->frozen);
    int i;
    double U = u*nnuc;
    double cost=0.0;
    double obval;  // observed mutations
    double exval;  // mutations expected under model
    for(i=0; i < BT_DIM; ++i) {
        BTLink *o = obs->tab[i];
        BTLink *e = expt->tab[i];
        while(o && e) {
            if(o->key < e->key) {
                // e link missing, so exval=0.
                obval = o->value;
                exval = 0.0;
                o = o->next;
            }else if(o->key > e->key) {
                // o link missing, so obval=0.
                obval = 0.0;
                exval = e->value * U;
                e = e->next;
            }else {
                assert(o->key == e->key);
                obval = o->value;
                exval = e->value * U;
                e = e->next;
                o = o->next;
            }
            cost += -obval*log(exval) + exval + lgamma(obval+1.0);
        }
        if(o) { // e link missing, so exval=0
            cost = HUGE_VAL;
            break;
        }
        obval=0.0;
        while(e) { // o link missing, so obval=0
            exval = e->value * U;
            cost += -obval*log(exval) + exval + lgamma(obval+1.0);
            e = e->next;
        }
    }
    assert(cost >= 0.0);
    return cost;
}
#endif

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
                // e->value is q=0. This case blows up unless p==0.
                p = o->value;
                if(p != 0.0) {
                    kl = HUGE_VAL;
#if 0
                    fprintf(stderr,"%s:%s:%d: missing expt for ",
                            __FILE__,__func__,__LINE__);
                    printBits(sizeof(o->key), &o->key, stderr);
                    //exit(EXIT_FAILURE);
#endif
                    o = o->next;
                } else {
                    o = o->next;
                    continue;
                }
            }else if(o->key > e->key) {
                // o->value is p=0. For this case, note that
                // p*log(p/q) is the log of (p^p / q^p), which -> 1 as
                // p->0.  The log of 1 is 0, so 0 is the contribution
                // to kl.
                e = e->next;
                continue;
            }else {
                assert(o->key == e->key);
                p = o->value;
                q = e->value;
                e = e->next;
                o = o->next;
                kl += p*log(p/q);
            }
        }
        while(o) { // e->value is q=0: fail unless p=0
            p = o->value;
            if(p != 0.0) {
                kl = HUGE_VAL;
#if 0
                fprintf(stderr,"%s:%s:%d: missing expt for ",
                        __FILE__,__func__,__LINE__);
                printBits(sizeof(o->key), &o->key, stderr);
                exit(EXIT_FAILURE);
#endif
                o = o->next;
            } else {
                o = o->next;
                continue;
            }
        }
        // Any remaining cases have p=0, so contribution to kl is 0.
    }
    return kl;
}

#ifdef TEST

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
    "twoN free   2Na=100\n"
    "twoN fixed  2Nb=123\n"
    "twoN free   2Nc=213.4\n"
    "twoN fixed  2Nbb=32.1\n"
    "twoN free   2Nab=222\n"
    "twoN fixed  2Nabc=1.2e2\n"
    "mixFrac free Mc=0.02\n"
    "segment a   t=T0     twoN=2Na    samples=1\n"
    "segment b   t=T0     twoN=2Nb    samples=1\n"
    "segment c   t=Tc     twoN=2Nc    samples=1\n"
    "segment bb  t=Tc     twoN=2Nbb\n"
    "segment ab  t=Tab    twoN=2Nab\n"
    "segment abc t=Tabc   twoN=2Nabc\n"
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

    if(verbose)
        BranchTab_print(bt);
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

    BranchTab *bt2 = BranchTab_dup(bt);
    assert(BranchTab_equals(bt, bt2));

    if(verbose)
        BranchTab_print(bt);
    BranchTab_free(bt);
    GPTree_free(g);
	unitTstResult("BranchTab", "OK");
    unlink(tstFname);
    unlink(tstPatProbFname);
}
#endif
