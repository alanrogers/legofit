#include "gptree.h"
#include "hashtab.h"
#include "misc.h"
#include "parse.h"
#include "sampndx.h"
#include "parstore.h"
#include "tokenizer.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Consider the following tree of populations:
//
//      a-------|
//              |ab--|
//      b--|bb--|    |
//         |         |abc--
//         |c--------|
//
//  t = 0  1    3    5.5     inf
//
// In this tree, a, b, c, bb, ab, and abc represent segments of the
// population tree.  The input file begins with a series of "segment"
// statements that define each of the segments in the tree. The
// segment statements also provide the time (backwards from the
// present) at which the segment starts, and the size, twoN, of the
// population. Optionally, it also provides the number of haploid
// samples observed in this segment of the tree.
//
// The statements that follow the segment statements describe how the
// segments are connected. The "mix" statement is used when a segment
// originates as a mixture of two ancestral segments. The "derive"
// statement is used when a segment derives from a single ancestral
// segment. 
//
// No segment can have more than two "parents" or more than two
// "children".
//
// Here is input that would generate the tree above:
//
// segment a   t=0     twoN=100    samples=1
// segment b   t=0     twoN=123    samples=2
// segment c   t=1     twoN=213.4  samples=1
// segment bb  t=1     twoN=32.1
// segment ab  t=3     twoN=222
// segment abc t=5.5e0 twoN=1.2e2
// mix    b  from bb + 0.1 c
// derive a  from ab
// derive bb from ab
// derive ab from abc
// derive c  from abc

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

#define ILLEGAL_INPUT(x) do{                                    \
        fprintf(stderr,"%s:%s:%d: Illegal input: \"%s\"\n",     \
            __FILE__,__func__,__LINE__, (x));                   \
        exit(EXIT_FAILURE);                                     \
    }while(0)

int         getDbl(double *x, Tokenizer * tkz, int i);
int         getULong(unsigned long *x, Tokenizer * tkz, int i);
void		parseParam(Tokenizer *tkz, ParKeyVal **pkv, ParStore *fixed,
					   ParStore *var);
void        parseSegment(Tokenizer *tkz, HashTab *poptbl, ParKeyVal **pkv,
						 SampNdx *sndx);
void        parseDerive(Tokenizer *tkz, HashTab *poptbl); 
void        parseMix(Tokenizer *tkz, HashTab *poptbl, ParKeyVal **pkv);

int getDbl(double *x, Tokenizer * tkz, int i) {
    char       *end = NULL;
	const char *tok = Tokenizer_token(tkz, i);

    *x = strtod(tok, &end);
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

// param fixed|varies name=100
void		parseParam(Tokenizer *tkz, ParKeyVal **pkv, ParStore *fixed,
					   ParStore *varies) {
	int curr=1, ntokens = Tokenizer_ntokens(tkz);
	int isfixed=0;

	// Read type of parameter: "fixed" or "varies"
	{
		char *tok;
		CHECK_INDEX(curr, ntokens);
		tok = Tokenizer_token(tkz, curr++);
		if(0 == strcmp("fixed", tok))
			fixed = 1;
		else if(0 == strcmp("varies", tok))
			fixed = 0;
		else {
			fprintf(stderr, "%s:%s:%d: got %s when expecting"
					" \"fixed\" or \"varies\".\n",
					__FILE__,__func__,__LINE__, tok);
			Tokenizer_print(tkz, stderr);
			exit(EXIT_FAILURE);
		}
	}

	// Read parameter name
	char *name;
	CHECK_INDEX(curr, ntokens);
	name = Tokenizer_token(tkz, curr++);
	assert(name != NULL);
	if( !isalpha(*name) )
		eprintf("%s:%s:%d: \"%s\" is not a legal parameter name.\n",
				__FILE__,__func__,__LINE__, name);

	// Read parameter value
    CHECK_INDEX(curr, ntokens);
	double value;
    if(getDbl(&value, tkz, curr)) {
		fflush(stdout);
        fprintf(stderr, 
                "%s:%s:%d:Can't parse \"%s\" as a double.\n",
				__FILE__,__func__,__LINE__,
                Tokenizer_token(tkz, curr));
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    }

	// Allocate and initialize parameter in ParStore
	double ptr;
	if(isfixed)
		ptr = ParStore_addPar(fixed, value);
	else
		ptr = ParStore_addPar(varies, value);

	// Add parameter name and pointer to linked list
	*pkv = ParKeyVal_add(*pkv, name, ptr);
}

// segment a   t=0     twoN=100    samples=1
void parseSegment(Tokenizer *tkz, HashTab *poptbl, SampNdx *sndx,
				  ParStore *fixed, ParStore *var, Bounds *bnd) {
    char *popName;
    double t, twoN;
    unsigned long nsamples=0;
    int curr=1,  ntokens = Tokenizer_ntokens(tkz);

    // Read name of segment
    CHECK_INDEX(curr, ntokens);
    popName = Tokenizer_token(tkz, curr++);

    // Read t
    CHECK_INDEX(curr, ntokens);
    if(0 != strcmp("t", Tokenizer_token(tkz, curr++))) {
        fprintf(stderr, "Got %s when expecting t on input:\n",
                Tokenizer_token(tkz, curr - 1));
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    }
    CHECK_INDEX(curr, ntokens);
    if(getDbl(&t, tkz, curr)) {
		fflush(stdout);
        fprintf(stderr, 
                "Can't parse \"%s\" as a double. Expecting value of t\n",
                Tokenizer_token(tkz, curr));
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    } else {
        if(t < 0.0) {
            fprintf(stderr, "Bad value of t: \"%s\" in input line:\n",
                    Tokenizer_token(tkz, curr));
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        ++curr;
    }
    
    // Read twoN
	if(curr >= ntokens) {
		fprintf(stderr, "curr=%d >= ntokens=%d\n", curr, ntokens);
		Tokenizer_print(tkz, stderr);
	}
    CHECK_INDEX(curr, ntokens);
    if(0 != strcmp("twoN", Tokenizer_token(tkz, curr++))) {
        fprintf(stderr, "Got %s when expecting twoN on input:\n",
                Tokenizer_token(tkz, curr - 1));
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    }
    CHECK_INDEX(curr, ntokens);
	int vary_twoN;
    if(getDbl(&twoN, &vary_twoN, tkz, curr)) {
        fprintf(stderr,
                "Can't parse \"%s\" as a double. Expecting value of twoN\n",
                Tokenizer_token(tkz, curr));
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    } else {
        if(twoN <= 0.0) {
            fprintf(stderr, "Bad value of twoN: \"%s\" in input line:\n",
                    Tokenizer_token(tkz, curr));
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        ++curr;
    }

    // Read (optional) number of samples
    if(curr < ntokens) {
        if(0 != strcmp("samples", Tokenizer_token(tkz, curr++)))
            eprintf("%s:%s:%d: got %s when expecting \"samples\"\n",
                     __FILE__,__func__,__LINE__,
                     Tokenizer_token(tkz, curr - 1));
        CHECK_INDEX(curr, ntokens);
        if(getULong(&nsamples, tkz, curr)) 
            eprintf("%s:%s:%d: Can't parse \"%s\" as an unsigned int."
                    " Expecting value of \"samples\"\n",
                     __FILE__,__func__,__LINE__,
                     Tokenizer_token(tkz, curr - 1));
        else {
            if(nsamples > MAXSAMP) 
                eprintf("%s:%s:%d: %lu samples is too many: max is %d:\n",
                         __FILE__,__func__,__LINE__, nsamples, MAXSAMP);
            ++curr;
        }
    }

    if(curr < ntokens)
        eprintf("%s:%s:%d: extra token \"%s\" at end of line\n",
                 __FILE__,__func__,__LINE__, Tokenizer_token(tkz,curr));

    assert(strlen(popName) > 0);
    El *e = HashTab_get(poptbl, popName);
    PopNode *thisNode = (PopNode *) El_get(e);
    if(thisNode == NULL) {
		// Allocate storage for parameters within ParStore.
		double *twoNptr, *startPtr, *endPtr, *mPtr;

		ParStore *parstore = (vary_twoN ? var : fixed);
		twoNptr = ParStore_addPar(parstore, twoN, bnd->lo_twoN, bnd->hi_twoN);

		parstore = (vary_t ? var : fixed);
		startPtr = ParStore_addPar(var, t, bnd->lo_t, bnd->hi_t);

		endPtr = ParStore_addPar(fixed, HUGE_VAL, bnd->lo_t, HUGE_VAL);

		// Problem: mPtr is allocated here, but I don't know yet
		// whether it is fixed or variable.
		mPtr = ParStore_addPar(fixed, 0.0, 0.0, 1.0);

		// Create new node, with pointers to the newly-allocated parameters
        thisNode = PopNode_new(twoNptr, startPtr, endPtr, mPtr);
        El_set(e, thisNode);
    }else
        eprintf("%s:%s:%d: duplicate \"segment %s\"\n",
                 __FILE__,__func__,__LINE__, popName);
    SampNdx_addSamples(sndx, nsamples, thisNode, popName);
}

// derive a from ab
void parseDerive(Tokenizer *tkz, HashTab *poptbl) {
    char *childName, *parName;
    int curr=1,  ntokens = Tokenizer_ntokens(tkz);

    // Read name of child
    CHECK_INDEX(curr, ntokens);
    childName = Tokenizer_token(tkz, curr++);

    // Read "from"
    CHECK_INDEX(curr, ntokens);
    if(0 != strcmp("from", Tokenizer_token(tkz, curr++))) {
        fprintf(stderr, "Got %s when expecting \"from\" on input:\n",
                Tokenizer_token(tkz, curr - 1));
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    }
    
    // Read name of parent
    CHECK_INDEX(curr, ntokens);
    parName = Tokenizer_token(tkz, curr++);

    if(curr < ntokens)
        eprintf("%s:%s:%d: extra token \"%s\" at end of line\n",
                 __FILE__,__func__,__LINE__, Tokenizer_token(tkz,curr));

    assert(strlen(childName) > 0);
    El *childEl = HashTab_get(poptbl, childName);
    PopNode *childNode = (PopNode *) El_get(childEl);
    if(childNode == NULL) {
        eprintf("%s:%s:%d: child segment \"%s\" undefined\n",
                 __FILE__,__func__,__LINE__, childName);
    }

    assert(strlen(parName) > 0);
    El *parEl = HashTab_get(poptbl, parName);
    PopNode *parNode = (PopNode *) El_get(parEl);
    if(parNode == NULL) {
        eprintf("%s:%s:%d: parent segment \"%s\" undefined\n",
                 __FILE__,__func__,__LINE__, parName);
    }
    PopNode_addChild(parNode, childNode);
}

// mix    b  from bb + 0.1 c
void parseMix(Tokenizer *tkz, HashTab *poptbl) {
    char *childName, *parName[2];
    double m;
    int curr=1,  ntokens = Tokenizer_ntokens(tkz);

    // Read name of child
    CHECK_INDEX(curr, ntokens);
    childName = Tokenizer_token(tkz, curr++);

    // Read word "from"
    CHECK_INDEX(curr, ntokens);
    if(0 != strcmp("from", Tokenizer_token(tkz, curr++)))
        eprintf("%s:%s:%d: got %s when expecting \"from\" on input:\n",
                 __FILE__,__func__,__LINE__, Tokenizer_token(tkz, curr - 1));

    // Read name of parent0
    CHECK_INDEX(curr, ntokens);
    parName[0] = Tokenizer_token(tkz, curr++);

    // Read symbol "+"
    CHECK_INDEX(curr, ntokens);
    if(0 != strcmp("+", Tokenizer_token(tkz, curr++)))
        eprintf("%s:%s:%d: got %s when expecting \"+\" on input:\n",
                 __FILE__,__func__,__LINE__, Tokenizer_token(tkz, curr - 1));

    // Read mixture fraction
    CHECK_INDEX(curr, ntokens);
    if(getDbl(&m, NULL, tkz, curr++) || m < 0.0 || m > 1.0 ) {
        eprintf("%s:%s:%d: bad mixture fraction \"%s\"->%0.20lf\n",
                 __FILE__,__func__,__LINE__,
                 Tokenizer_token(tkz, curr-1), m);
    }
    // Read name of parent1
    CHECK_INDEX(curr, ntokens);
    parName[1] = Tokenizer_token(tkz, curr++);

    if(curr < ntokens)
        eprintf("%s:%s:%d: extra token \"%s\" at end of line\n",
                 __FILE__,__func__,__LINE__, Tokenizer_token(tkz,curr));

    assert(m <= 1.0);
    assert(m >= 0.0);

    assert(strlen(childName) > 0);
    El *childEl = HashTab_get(poptbl, childName);
    PopNode *childNode = (PopNode *) El_get(childEl);
    if(childNode == NULL) {
        eprintf("%s:%s:%d: child segment \"%s\" undefined\n",
                 __FILE__,__func__,__LINE__, childName);
    }

    assert(strlen(parName[0]) > 0);
    El *parEl0 = HashTab_get(poptbl, parName[0]);
    PopNode *parNode0 = (PopNode *) El_get(parEl0);
    if(parNode0 == NULL)
        eprintf("%s:%s:%d: parent segment \"%s\" undefined\n",
                 __FILE__,__func__,__LINE__, parName[0]);

    assert(strlen(parName[1]) > 0);
    El *parEl1 = HashTab_get(poptbl, parName[1]);
    PopNode *parNode1 = (PopNode *) El_get(parEl1);
    if(parNode1 == NULL)
        eprintf("%s:%s:%d: parent segment \"%s\" undefined\n",
                 __FILE__,__func__,__LINE__, parName[1]);

    PopNode_mix(childNode, m, parNode1, parNode0);
}

PopNode    *mktree(FILE * fp, HashTab * poptbl, SampNdx *sndx, ParStore *fixed,
				   ParStore *var, Bounds *bnd) {
    int         ntokens;
    char        buff[500];
    Tokenizer  *tkz = Tokenizer_new(50);
	ParKeyVal    *pkv = NULL;  // table of pointers to parameters

    while(1) {
        if(fgets(buff, sizeof(buff), fp) == NULL)
            break;

        if(!strchr(buff, '\n') && !feof(fp)) {
            fprintf(stderr, "ERR@%s:%d: input buffer overflow."
                    " buff size: %zu\n", __FILE__, __LINE__, sizeof(buff));
            exit(EXIT_FAILURE);
        }

        // strip trailing comments
        char *comment = strchr(buff, '#');
        if(comment)
            *comment = '\0';

        Tokenizer_split(tkz, buff, " \t="); // tokenize
        ntokens = Tokenizer_strip(tkz, " \t=\n");
        if(ntokens == 0)
            continue;

        char *tok = Tokenizer_token(tkz, 0);
		if(0 == strcmp(tok, "param"))
			parseParam(tkz, &pkv, fixed, var, bnd);
        else if(0 == strcmp(tok, "segment"))
            parseSegment(tkz, poptbl, &pkv, sndx);
        else if(0 == strcmp(tok, "mix"))
            parseMix(tkz, poptbl, &pkv);
        else if(0 == strcmp(tok, "derive"))
            parseDerive(tkz, poptbl);
        else
            ILLEGAL_INPUT(tok);
    }

    // Make sure the tree of populations has a single root. This
    // code iterates through all the nodes in the HashTab, and
    // searches from each node back to the root. If all is well,
    // these searches all find the same root. Otherwise, it aborts
    // with an error.
    PopNode    *root = NULL;
    {
        HashTabSeq *popseq = HashTabSeq_new(poptbl);
        CHECKMEM(popseq);

        El         *el = HashTabSeq_next(popseq);
        while(el != NULL) {
            PopNode    *node = El_get(el);
            assert(node != NULL);
#if 0
            El_printShallow(el);
            PopNode_printShallow(node, stdout);
#endif

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
            el = HashTabSeq_next(popseq);
        }
        HashTabSeq_free(popseq);
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

//      a-------|
//              |ab--|
//      b--|bb--|    |
//         |         |abc--
//         |c--------|
//
//  t = 0  1    3    5.5     inf
const char *tstInput =
    " # this is a comment\n"
    "segment a   t=0     twoN=100    samples=1\n"
    "segment b   t=0     twoN=123    samples=2\n"
    "segment c   t=1     twoN=213.4  samples=1\n"
    "segment bb  t=1     twoN=32.1 # another comment\n"
    "segment ab  t=3     twoN=222\n"
    "segment abc t=5.5e0 twoN=1.2e2\n"
    "mix b from bb + 0.1 c\n"
    "derive a from ab\n"
    "derive bb from ab\n"
    "derive ab from abc\n"
    "derive c from abc\n";

int main(int argc, char **argv) {

    int verbose=0;

    if(argc > 1) {
        if(argc!=2 || 0!=strcmp(argv[1], "-v")) {
            fprintf(stderr,"usage: xparse [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

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

    HashTab *poptbl = HashTab_new();
    SampNdx sndx;
    SampNdx_init(&sndx);
	ParStore *fixed = ParStore_new();  // fixed parameters
	Bounds   bnd = {
		.lo_twoN = 0.0,
		.hi_twoN = 1e7,
		.lo_t = 0.0,
		.hi_t = HUGE_VAL
	};
    PopNode *root = mktree(fp, poptbl, &sndx, fixed, &bnd);

    if(verbose) {
        PopNode_print(stdout, root, 0);
        unsigned i;
        for(i=0; i < SampNdx_size(&sndx); ++i)
            printf("%2u %s\n", i, SampNdx_lbl(&sndx, i));
		printf("Used %d parameters in \"fixed\".\n", ParStore_nPar(fixed));
    }

    HashTab_free(poptbl);
	ParStore_free(fixed);
    return 0;
}
#endif
