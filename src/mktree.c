#include "gptree.h"
#include "hashtab.h"
#include "misc.h"
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
    double      start, end; // duration of interval in backwards time
    unsigned long nsamples;
};

#define CHECK_INDEX(ndx,n) do{                                  \
        if((ndx)>=(n)){                                         \
            fprintf(stderr,"%s:%s:%d: index out of bounds\n",   \
                    __FILE__,__func__,__LINE__);                \
            exit(EXIT_FAILURE);                                 \
        }                                                       \
    }while(0)

#define CHECK_NAME(s) do {                                      \
        if((s)==NULL) {                                         \
            fprintf(stderr, "%s:%s:%d: Null file name\n",       \
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

int         getDbl(double *x, Tokenizer * tkz, int i);
int         getULong(unsigned long *x, Tokenizer * tkz, int i);
PopData     parse(Tokenizer * tkz);
PopNode    *mktree(FILE * fp, HashTab * ht);

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
    double      m = 0, twoN = 0, start= 0, end=0;
    unsigned long nsamples = 0;
    int         curr = 0, ntokens = Tokenizer_ntokens(tkz);

    // Read name of current population
    CHECK_INDEX(curr, ntokens);
    thisPop = Tokenizer_token(tkz, curr++);

    // Read name of ancestral population
    CHECK_INDEX(curr, ntokens);
    if('<' != *Tokenizer_token(tkz, curr++)) {
        fprintf(stderr, "Bad input line. Got %s when expecting \"<\":\n",
                Tokenizer_token(tkz, curr - 1));
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    }
    CHECK_INDEX(curr, ntokens);
    ancestor = Tokenizer_token(tkz, curr++);

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
        char *s = Tokenizer_token(tkz, curr);
        char *eptr=NULL;
        start = strtod(s, &eptr);
        if(eptr==NULL || *eptr!=':') {
            fprintf(stderr,"%s:%s:%d: interval \"%s\" incorrectly formatted.\n",
                    __FILE__,__func__,__LINE__, s);
            fprintf(stderr,"          Should be like 123:452\n");
            exit(EXIT_FAILURE);
        }
        end = strtod(eptr+1, NULL);
        if(!(0.0 <= start && start<=end)) {
            fprintf(stderr,"%s:%s:%d: bad interval: \"%s\"\n",
                    __FILE__,__func__,__LINE__, s);
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

    CHECK_NAME(thisPop);
    snprintf(pd.thisPop, sizeof pd.thisPop, "%s", thisPop);
    CHECK_NAME(ancestor);
    snprintf(pd.ancestor, sizeof pd.ancestor, "%s", ancestor);
    if(m > 0.0) {
        pd.m = m;
        snprintf(pd.introgressor, sizeof pd.introgressor, "%s",
                 introgressor);
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
    char        buff[500];
    Tokenizer  *tkz = Tokenizer_new(50);
    PopData      pd;
    El     *thisEl, *ancestorEl, *introgressorEl;
    PopNode *thisNode, *ancestorNode, *introgressorNode;

    while(1) {
        if(fgets(buff, sizeof buff, fp) == NULL)
            break;

        if(!strchr(buff, '\n') && !feof(fp)) {
            fprintf(stderr, "ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, sizeof(buff));
            exit(EXIT_FAILURE);
        }
        Tokenizer_split(tkz, buff, " \t="); // tokenize
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
        }else
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
        }else
            PopNode_addChild(ancestorNode, thisNode);
    }

    // Make sure the tree of populations has a single root. This
    // code iterates through all the nodes in the HashTab, and
    // searches from each node back to the root. If all is well,
    // these searches all find the same root. Otherwise, it aborts
    // with an error.
    PopNode *root = NULL;
    {
        HashTabSeq *hts = HashTabSeq_new(ht);
        CHECKMEM(hts);

        El *el = HashTabSeq_next(hts);
        while(el != NULL) {
            PopNode *node = El_get(el);
            node = PopNode_root(node);
            if(root == NULL)
                root = node;
            else {
                if(root != node) {
                    fprintf(stderr,"%s:%s:%d: Pop tree has multiple roots.\n",
                            __FILE__,__func__,__LINE__);
                    exit(EXIT_FAILURE);
                }
            el = HashTabSeq_next(hts);
        }
        HashTabSeq_free(hts);
    }

    Tokenizer_free(tkz);

    return root;
}
