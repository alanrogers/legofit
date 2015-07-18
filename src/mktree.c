#include "gptree.h"
#include "hashtab.h"
#include "misc.h"
#include "mktree.h"
#include "tokenizer.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// One line of input for each segment of population tree.
// Each line is of form
//
//   x < y < 0.1 z 2N=1e5 t=1000:2000 samples=2
//
// which says that population x got 90% of its genes from y and 10%
// from z, has size 2N=1e5, during an interval from 1000 to 2000
// generations ago, and that we sampled 2 haploid genomes. If there
// was only a single parental population, the input line would look
// like this
//
//   x < y 2N=1e5 t=1000:2000 samples=2

#define POPNAMESIZE 30
#define MAXSAMPLES 10

struct PopData {
    char        thisPop[POPNAMESIZE];
    char        ancestor[POPNAMESIZE];
    char        introgressor[POPNAMESIZE];
    double      m, twoN;
    double      start, end;     // duration of interval in backwards time
    unsigned long nsamples;
};

#if 1
#  define CHECK_INDEX(ndx,n) do{                                \
        if((ndx)>=(n)){                                         \
            fprintf(stderr,"%s:%s:%d: index out of bounds\n",   \
                    __FILE__,__func__,__LINE__);                \
            exit(EXIT_FAILURE);                                 \
        }                                                       \
    }while(0)
#else
#  define CHECK_INDEX(ndx,n)
#endif

#if 1
#  define CHECK_NAME(s) do {                                    \
        if((s)==NULL) {                                         \
            fprintf(stderr, "%s:%s:%d: Null population name\n", \
                    __FILE__,__func__,__LINE__);                \
            exit(EXIT_FAILURE);                                 \
        }                                                       \
        if(strlen((s)) >= POPNAMESIZE) {                        \
            fprintf(stderr,"%s:%s:%d:"                          \
                    " Pop name \"%s\" is too long."             \
                    " Max=%d.\n",                               \
                    __FILE__,__func__,__LINE__,                 \
                    (s), POPNAMESIZE-1);                        \
            exit(EXIT_FAILURE);                                 \
        }                                                       \
    }while(0)
#else
#  define CHECK_NAME(s)
#endif

int         getDbl(double *x, Tokenizer * tkz, int i);
int         getULong(unsigned long *x, Tokenizer * tkz, int i);
PopData     parse(Tokenizer * tkz);
void addSpaces(size_t nTo, char to[nTo], size_t nFro, char from[nFro]);

void addSpaces(size_t nTo, char to[nTo], size_t nFro, char from[nFro]) {
    size_t i, j;

    i=j=0;
    while(j<nFro &&  from[j]!='\0') {
        if(from[j] == '<') {
            if(i > nTo-4) {
                fprintf(stderr,"%s:%s:%d: output buffer overflow\n",
                        __FILE__,__func__,__LINE__);
                exit(EXIT_FAILURE);
            }
            to[i++] = ' ';
            to[i++] = '<';
            to[i++] = ' ';
        }else{
            if(i > nTo-2) {
                fprintf(stderr,"%s:%s:%d: output buffer overflow\n",
                        __FILE__,__func__,__LINE__);
                exit(EXIT_FAILURE);
            }
            to[i++] = from[j];
        }
        ++j;
    }
    assert(i < nTo);
    to[i] = '\0';
}

int getDbl(double *x, Tokenizer * tkz, int i) {
    char       *end = NULL;

    *x = strtod(Tokenizer_token(tkz, i), &end);
    if(end != NULL && *end == '\0')
        return 0;               // success
    return 1;                   // failure
}

int getULong(unsigned long *x, Tokenizer * tkz, int i) {
    char       *end = NULL;

    *x = strtoul(Tokenizer_token(tkz, i), &end, 10);
    if(end != NULL && *end == '\0')
        return 0;               // success
    return 1;                   // failure
}

PopData parse(Tokenizer * tkz) {
    PopData     pd;
    memset(&pd, 0, sizeof pd);

    char       *thisPop = NULL, *ancestor = NULL, *introgressor = NULL;
    double      m = 0, twoN = 0, start = 0, end = 0;
    unsigned long nsamples = 0;
    int         curr = 0, ntokens = Tokenizer_ntokens(tkz);

    // Read name of current population
    CHECK_INDEX(curr, ntokens);
    thisPop = Tokenizer_token(tkz, curr++);

    // Read name of optional ancestral population
    CHECK_INDEX(curr, ntokens);
    if('<' == *Tokenizer_token(tkz, curr)) {
        ancestor = Tokenizer_token(tkz, ++curr);
        ++curr;
    }
    CHECK_INDEX(curr, ntokens);

    // Read admixture fraction and name of optional introgressing population
    CHECK_INDEX(curr, ntokens);
    if('<' == *Tokenizer_token(tkz, curr)) {
        if(getDbl(&m, tkz, ++curr)) {
            fprintf(stderr,
                    "Can't parse \"%s\" as a double. Expecting admixture fraction.\n",
                    Tokenizer_token(tkz, curr));
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        } else {
            if(m <= 0.0 || m >= 1.0) {
                fprintf(stderr, "Bad probability, %s, in input line:\n",
                        Tokenizer_token(tkz, curr - 1));
                Tokenizer_print(tkz, stderr);
                exit(EXIT_FAILURE);
            }
        }
        introgressor = Tokenizer_token(tkz, ++curr);
        ++curr;
    }

    // Read 2N
    CHECK_INDEX(curr, ntokens);
    if(0 != strcmp("2N", Tokenizer_token(tkz, curr++))) {
        fprintf(stderr, "Got %s when expecting 2N on input:\n",
                Tokenizer_token(tkz, curr - 1));
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    }
    CHECK_INDEX(curr, ntokens);
    if(getDbl(&twoN, tkz, curr)) {
        fprintf(stderr,
                "Can't parse \"%s\" as a double. Expecting value of 2N\n",
                Tokenizer_token(tkz, curr));
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    } else {
        if(twoN <= 0.0) {
            fprintf(stderr, "Bad value of 2N: \"%s\" in input line:\n",
                    Tokenizer_token(tkz, curr));
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        ++curr;
    }

    // Read duration
    CHECK_INDEX(curr, ntokens);
    if(0 != strcmp("t", Tokenizer_token(tkz, curr++))) {
        fprintf(stderr, "Got %s when expecting \"t\" on input:\n",
                Tokenizer_token(tkz, curr - 1));
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    }
    CHECK_INDEX(curr, ntokens);
    {
        char       *s = Tokenizer_token(tkz, curr);
        char       *eptr = NULL;
        start = strtod(s, &eptr);
        if(eptr == NULL || *eptr != ':') {
            fprintf(stderr,
                    "%s:%s:%d: interval \"%s\" incorrectly formatted.\n",
                    __FILE__, __func__, __LINE__, s);
            fprintf(stderr, "          Should be like 123:452\n");
            exit(EXIT_FAILURE);
        }
        end = strtod(eptr + 1, NULL);
        if(!(0.0 <= start && start <= end)) {
            fprintf(stderr, "%s:%s:%d: bad interval: \"%s\"\n",
                    __FILE__, __func__, __LINE__, s);
            exit(EXIT_FAILURE);
        }
        ++curr;
    }

    // Read (optional) number of samples
    if(curr < ntokens) {
        CHECK_INDEX(curr, ntokens);
        if(0 != strcmp("samples", Tokenizer_token(tkz, curr++))) {
            fprintf(stderr, "Got %s when expecting \"samples\" on input:\n",
                    Tokenizer_token(tkz, curr - 1));
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        CHECK_INDEX(curr, ntokens);
        if(getULong(&nsamples, tkz, curr)) {
            fprintf(stderr, "Can't parse \"%s\" as an unsigned int."
                    " Expecting value of \"samples\"\n",
                    Tokenizer_token(tkz, curr));
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        } else {
            if(nsamples > MAXSAMPLES) {
                fprintf(stderr, "Bad value samples=%lu. Max is %d:\n",
                        nsamples, MAXSAMPLES);
                Tokenizer_print(tkz, stderr);
                exit(EXIT_FAILURE);
            }
            ++curr;
        }
    }

    if(curr < ntokens) {
        fprintf(stderr, "Too many tokens on input line:\n");
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    }

    if(ancestor==NULL && !isinf(end)) {
        fprintf(stderr, "%s:%s:%d: ancestor is missing but interval is finite\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }
    if(ancestor!=NULL && isinf(end)) {
        fprintf(stderr, "%s:%s:%d: ancestor is present but interval is infinite\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }

    CHECK_NAME(thisPop);
    snprintf(pd.thisPop, sizeof pd.thisPop, "%s", thisPop);
    if(ancestor) {
        CHECK_NAME(ancestor);
        snprintf(pd.ancestor, sizeof pd.ancestor, "%s", ancestor);
    }
    if(m > 0.0) {
        pd.m = m;
        snprintf(pd.introgressor, sizeof pd.introgressor, "%s", introgressor);
    }
    pd.twoN = twoN;
    pd.start = start;
    pd.end = end;
    if(nsamples > 0)
        pd.nsamples = nsamples;

    return pd;
}

PopNode    *mktree(FILE * fp, HashTab * ht) {
    int         ntokens;
    char        buff0[500];
    char        buff1[500];
    Tokenizer  *tkz = Tokenizer_new(50);
    PopData     pd;
    El         *thisEl, *ancestorEl, *introgressorEl;
    PopNode    *thisNode, *ancestorNode, *introgressorNode;

    while(1) {
        if(fgets(buff0, sizeof buff0, fp) == NULL)
            break;

        if(!strchr(buff0, '\n') && !feof(fp)) {
            fprintf(stderr, "ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, sizeof(buff0));
            exit(EXIT_FAILURE);
        }
        addSpaces(sizeof buff1, buff1, sizeof buff0, buff0);
        Tokenizer_split(tkz, buff1, " \t="); // tokenize
        ntokens = Tokenizer_strip(tkz, " \t=\n");
        if(ntokens == 0)
            continue;

        pd = parse(tkz);

        thisEl = HashTab_get(ht, pd.thisPop);
        thisNode = (PopNode *) El_get(thisEl);
        if(thisNode == NULL) {
            thisNode = PopNode_new(pd.twoN, pd.start, pd.end);
            CHECKMEM(thisNode);
            El_set(thisEl, thisNode);
        } else
            PopNode_set(thisNode, pd.twoN, pd.start, pd.end);

        ancestorEl = HashTab_get(ht, pd.ancestor);
        ancestorNode = (PopNode *) El_get(ancestorEl);
        if(ancestorNode == NULL) {
            ancestorNode = PopNode_new(-1.0, -1.0, -1.0);
            CHECKMEM(ancestorNode);
            El_set(ancestorEl, ancestorNode);
        }

        if(pd.m > 0) {
            introgressorEl = HashTab_get(ht, pd.introgressor);
            introgressorNode = (PopNode *) El_get(introgressorEl);
            if(introgressorNode == NULL) {
                introgressorNode = PopNode_new(-1.0, -1.0, -1.0);
                CHECKMEM(introgressorNode);
                El_set(introgressorEl, introgressorNode);
            }

            PopNode_mix(thisNode, pd.m, introgressorNode, ancestorNode);
        } else
            PopNode_addChild(ancestorNode, thisNode);
    }

    // Make sure the tree of populations has a single root. This
    // code iterates through all the nodes in the HashTab, and
    // searches from each node back to the root. If all is well,
    // these searches all find the same root. Otherwise, it aborts
    // with an error.
    PopNode    *root = NULL;
    {
        HashTabSeq *hts = HashTabSeq_new(ht);
        CHECKMEM(hts);

        El         *el = HashTabSeq_next(hts);
        while(el != NULL) {
            PopNode    *node = El_get(el);
            assert(node != NULL);
            El_print(el);
            PopNode_printShallow(node, stdout);
#if 0
            PopNode_sanityFromLeaf(node, __FILE__, __LINE__);
            node = PopNode_root(node);
            if(root == NULL)
                root = node;
            else {
                if(root != node) {
                    fprintf(stderr,
                            "%s:%s:%d: Pop tree has multiple roots.\n",
                            __FILE__, __func__, __LINE__);
                    exit(EXIT_FAILURE);
                }
            }
#endif
            el = HashTabSeq_next(hts);
        }
        HashTabSeq_free(hts);
    }
    Tokenizer_free(tkz);
    return root;
}

#ifdef TEST

#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

PopData tryparse(const char *s, int verbose);

//  a
//    > ab---  abc
//  b        /
//      /----
//  c -
//
const char *tstInput = 
    "a  < ab           2N=1e5 t=0:10 samples=2\n"
    "b  < bb < 0.25 bc 2N=123 t=0:5  samples=1\n"
    "bb < ab           2N=44  t=5:10          \n"
    "c  < bc           2N=234 t=0:5  samples=1\n"
    "bc < abc          2N=1   t=5:10          \n"
    "ab < abc          2N=888 t=10:20         \n"
    "abc               2N=2   t=10:inf        \n";

PopData tryparse(const char *s, int verbose) {
    Tokenizer *tkz = Tokenizer_new(50);
    const int bsize = 100;
    char buff0[bsize];
    char buff1[bsize];
    int ntokens;
    PopData pd;

    snprintf(buff0, sizeof buff0, "%s", s);
    if(verbose)
        printf("buff0: %s\n", buff0);
    addSpaces(sizeof buff1, buff1, sizeof buff0, buff0);
    if(verbose)
        printf("buff1: %s\n", buff1);
    Tokenizer_split(tkz, buff1, " \t=");
    ntokens = Tokenizer_strip(tkz, " \t=\n");
    if(verbose)
        printf("ntokens=%d\n", ntokens);
    pd = parse(tkz);
    if(verbose) {
        Tokenizer_print(tkz, stdout);
        printf("ego=%s parent=%s introg=%s m=%lf twoN=%lf t=%lf:%lf nsamples=%lu\n",
               pd.thisPop, pd.ancestor, pd.introgressor,
               pd.m, pd.twoN, pd.start, pd.end, pd.nsamples);
        putchar('\n');
    }

    Tokenizer_free(tkz);
    return pd;
}

int main(int argc, char **argv) {

    int verbose=0;
    PopData pd;

    if(argc > 1) {
        if(argc!=2 || 0!=strcmp(argv[1], "-v")) {
            fprintf(stderr,"usage: xmktree [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    pd = tryparse("x < yy 2N=12e3 t=100:200", verbose);
    assert(0 == strcmp(pd.thisPop, "x"));
    assert(0 == strcmp(pd.ancestor, "yy"));
    assert(0 == strlen(pd.introgressor));
    assert(0.0 == pd.m);
    assert(12e3 == pd.twoN);
    assert(100 == pd.start);
    assert(200 == pd.end);
    assert(0L == pd.nsamples);

    pd = tryparse("x < y 2N=12e3 t=100:200 samples=3", verbose);
    assert(0 == strcmp(pd.thisPop, "x"));
    assert(0 == strcmp(pd.ancestor, "y"));
    assert(0 == strlen(pd.introgressor));
    assert(0.0 == pd.m);
    assert(12e3 == pd.twoN);
    assert(100 == pd.start);
    assert(200 == pd.end);
    assert(3 == pd.nsamples);

    pd = tryparse("x < y < 0.25 z 2N=12e3 t=100:200 samples=3", verbose);
    assert(0 == strcmp(pd.thisPop, "x"));
    assert(0 == strcmp(pd.ancestor, "y"));
    assert(0 == strcmp(pd.introgressor, "z"));
    assert(0.25 == pd.m);
    assert(12e3 == pd.twoN);
    assert(100 == pd.start);
    assert(200 == pd.end);
    assert(3 == pd.nsamples);

    pd = tryparse("x<y<0.25 z 2N = 12e3 t=100:200 samples = 4", verbose);
    assert(0 == strcmp(pd.thisPop, "x"));
    assert(0 == strcmp(pd.ancestor, "y"));
    assert(0 == strcmp(pd.introgressor, "z"));
    assert(0.25 == pd.m);
    assert(12e3 == pd.twoN);
    assert(100 == pd.start);
    assert(200 == pd.end);
    assert(4 == pd.nsamples);

    const char *tstFname = "mktree-tmp.lgo";
    FILE       *fp = fopen(tstFname, "w");
    fputs(tstInput, fp);
    fclose(fp);
    fp = fopen(tstFname, "r");
    if(fp == NULL) {
        fprintf(stderr, "%s:%d: Can't open file \"%s\"\n",
                __FILE__, __LINE__, tstFname);
        exit(1);
    }

    HashTab *ht = HashTab_new();
    PopNode *root = mktree(fp, ht);

    HashTab_free(ht);
    return 0;
}
#endif
