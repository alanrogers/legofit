#include "exopar.h"
#include "popnodetab.h"
#include "lblndx.h"
#include "misc.h"
#include "parse.h"
#include "parstore.h"
#include "popnode.h"
#include "tokenizer.h"
#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

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
// time fixed  T0=0
// time free Tc=1
// time free Tab=3
// time free Tabc=5.5
// twoN free 2Na=100
// twoN fixed  2Nb=123
// twoN free 2Nc=213.4
// twoN fixed  2Nbb=32.1
// twoN free 2Nab=222
// twoN fixed  2Nabc=1.2e2
// mixFrac free Mc=0.8
// segment a   t=T0     twoN=2Na    samples=1
// segment b   t=T0     twoN=2Nb    samples=2
// segment c   t=Tc     twoN=2Nc    samples=1
// segment bb  t=Tc     twoN=2Nbb
// segment ab  t=Tab    twoN=2Nab
// segment abc t=Tabc   twoN=2Nabc
// mix    b  from bb + Mc * c
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

enum ParamType { TwoN, Time, MixFrac };
enum ParamStatus { Free, Fixed, Gaussian };

int         getDbl(double *x, Tokenizer * tkz, int i);
int         getULong(unsigned long *x, Tokenizer * tkz, int i);
void		parseParam(Tokenizer *tkz, enum ParamType type,
					   ParStore *parstore, Bounds *bnd,
                       ExoPar *exopar, bool isTimeParam);
void        parseSegment(Tokenizer *tkz, PopNodeTab *poptbl, SampNdx *sndx,
						 LblNdx *lndx, ParStore *parstore,
                         NodeStore *ns);
void        parseDerive(Tokenizer *tkz, PopNodeTab *poptbl);
void        parseMix(Tokenizer *tkz, PopNodeTab *poptbl, ParStore *parstore);

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

// twoN    fixed|free name=100
// time    fixed|free name=100
// mixFrac fixed|free name=100
void		parseParam(Tokenizer *tkz, enum ParamType type,
					   ParStore *parstore, Bounds *bnd,
                       ExoPar *exopar, bool isTimeParam) {
	int curr=1, ntokens = Tokenizer_ntokens(tkz);
    enum ParamStatus pstat = Free;

	// Read type of parameter: "fixed" or "free"
	{
		char *tok;
		CHECK_INDEX(curr, ntokens);
		tok = Tokenizer_token(tkz, curr++);
		if(0 == strcmp("fixed", tok))
			pstat = Fixed;
		else if(0 == strcmp("free", tok))
			pstat = Free;
        else if(0 == strcmp("gaussian", tok))
            pstat = Gaussian;
		else {
			fprintf(stderr, "%s:%s:%d: got %s when expecting"
					" \"fixed\" or \"free\".\n",
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
	if( !isalnum(*name) )
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
    ++curr;

    double sd = 0.0;
    if(pstat == Gaussian) {

        // Read string "sd"
        char *sdstr;
        CHECK_INDEX(curr, ntokens);
        sdstr = Tokenizer_token(tkz, curr++);
        assert(sdstr != NULL);
        if( 0 != strcmp(sdstr, "sd") )
            eprintf("%s:%s:%d: \"%s\" should be \"sd\".\n",
                    __FILE__,__func__,__LINE__, sdstr);

        // Read sd value
        CHECK_INDEX(curr, ntokens);
        if(getDbl(&sd, tkz, curr)) {
            fflush(stdout);
            fprintf(stderr,
                    "%s:%s:%d:Can't parse \"%s\" as a double.\n",
                    __FILE__,__func__,__LINE__,
                    Tokenizer_token(tkz, curr));
            Tokenizer_print(tkz, stderr);
            exit(EXIT_FAILURE);
        }
        if(sd <= 0.0) {
            fprintf(stderr,"%s:%s:%d: Warning: sd <= 0.0\n",
                    __FILE__,__func__,__LINE__);
            Tokenizer_print(tkz, stderr);
        }
        ++curr;
    }

	// Allocate and initialize parameter in ParStore
    switch(pstat) {
    case Fixed:
		ParStore_addFixedPar(parstore, value, name);
        break;
    case Gaussian:
		ParStore_addFixedPar(parstore, value, name);
        {
            bool isfree;
            double *ptr = ParStore_findPtr(parstore, &isfree, name);
            if(ptr == NULL) {
                fprintf(stderr,"%s:%d: Parameter \"%s\" is undefined.\n",
                        __FILE__,__LINE__, name);
                Tokenizer_print(tkz, stderr);
                exit(EXIT_FAILURE);
            }
            assert(!isfree);
            ExoPar_add(exopar, ptr, value, sd);
        }
        break;
    case Free:
        // Set bounds, based on type of parameter
        {
            double lo, hi;
            switch(type) {
            case TwoN:
                lo = bnd->lo_twoN;
                hi = bnd->hi_twoN;
                break;
            case Time:
                lo = bnd->lo_t;
                hi = bnd->hi_t;
                break;
            case MixFrac:
                lo = 0.0;
                hi = 1.0;
                break;
            default:
                lo = hi = 0.0;
                eprintf("%s:%s:%d: This shouldn't happen\n",
                        __FILE__,__func__,__LINE__);
            }
            ParStore_addFreePar(parstore, value, lo, hi, name, isTimeParam);
        }
        break;
    default:
        eprintf("%s:%s:%d: This shouldn't happen\n",
                __FILE__,__func__,__LINE__);
    }
}

// segment a   t=0     twoN=100    samples=1
void parseSegment(Tokenizer *tkz, PopNodeTab *poptbl, SampNdx *sndx,
				  LblNdx *lndx, ParStore *parstore, NodeStore *ns) {
    char *popName, *tok;
    double *tPtr, *twoNptr;
	bool tfree, twoNfree;
    unsigned long nsamples=0;
    int curr=1,  ntokens = Tokenizer_ntokens(tkz);

    // Read name of segment
    CHECK_INDEX(curr, ntokens);
    popName = Tokenizer_token(tkz, curr++);

    // Read t
    CHECK_INDEX(curr, ntokens);
    if(0 != strcmp("t", Tokenizer_token(tkz, curr++))) {
        fprintf(stderr, "Got %s when expecting \"t\" on input:\n",
                Tokenizer_token(tkz, curr - 1));
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    }
    CHECK_INDEX(curr, ntokens);
    tok = Tokenizer_token(tkz, curr++);
    tPtr = ParStore_findPtr(parstore, &tfree, tok);
	if(NULL == tPtr) {
		fprintf(stderr,"%s:%s:%d: Parameter \"%s\" is undefined\n",
				__FILE__,__func__,__LINE__,tok);
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    }

    // Read twoN
	if(curr >= ntokens) {
		fprintf(stderr, "curr=%d >= ntokens=%d\n", curr, ntokens);
		Tokenizer_print(tkz, stderr);
	}
    CHECK_INDEX(curr, ntokens);
    if(0 != strcmp("twoN", Tokenizer_token(tkz, curr++))) {
        fprintf(stderr, "Got %s when expecting \"twoN\" on input:\n",
                Tokenizer_token(tkz, curr - 1));
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    }
    CHECK_INDEX(curr, ntokens);
    tok = Tokenizer_token(tkz, curr++);
    twoNptr = ParStore_findPtr(parstore, &twoNfree, tok);
	if(NULL == twoNptr) {
		fprintf(stderr,"%s:%s:%dParameter \"%s\" is undefined\n",
				__FILE__,__func__,__LINE__, tok);
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
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
    PopNode *thisNode = PopNode_new(twoNptr, twoNfree, tPtr, tfree, ns);
    if(0 != PopNodeTab_insert(poptbl, popName, thisNode))
        eprintf("%s:%s:%d: duplicate \"segment %s\"\n",
                 __FILE__,__func__,__LINE__, popName);
    LblNdx_addSamples(lndx, nsamples, popName);
    SampNdx_addSamples(sndx, nsamples, thisNode);
}

// derive a from ab
void parseDerive(Tokenizer *tkz, PopNodeTab *poptbl) {
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
    PopNode *childNode = PopNodeTab_get(poptbl, childName);
    if(childNode == NULL) {
        eprintf("%s:%s:%d: child segment \"%s\" undefined\n",
                 __FILE__,__func__,__LINE__, childName);
    }

    assert(strlen(parName) > 0);
    PopNode *parNode = PopNodeTab_get(poptbl, parName);
    if(parNode == NULL) {
        eprintf("%s:%s:%d: parent segment \"%s\" undefined\n",
                 __FILE__,__func__,__LINE__, parName);
    }
    PopNode_addChild(parNode, childNode);
}

// mix    b  from bb + 0.1 c
void parseMix(Tokenizer *tkz, PopNodeTab *poptbl, ParStore *parstore) {
    char *childName, *parName[2], *tok;
    double *mPtr;
	bool mfree;
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

    // Read mixture fraction, mPtr
    CHECK_INDEX(curr, ntokens);
    CHECK_INDEX(curr, ntokens);
    tok = Tokenizer_token(tkz, curr++);
    mPtr = ParStore_findPtr(parstore, &mfree, tok);
	if(NULL == mPtr) {
		fprintf(stderr,"%s:%s:%d: Parameter \"%s\" is undefined\n",
				__FILE__,__func__,__LINE__, tok);
        Tokenizer_print(tkz, stderr);
        exit(EXIT_FAILURE);
    }

    // Read symbol "*"
    CHECK_INDEX(curr, ntokens);
    if(0 != strcmp("*", Tokenizer_token(tkz, curr++)))
        eprintf("%s:%s:%d: got %s when expecting \"*\" on input:\n",
                 __FILE__,__func__,__LINE__, Tokenizer_token(tkz, curr - 1));

    // Read name of parent1
    CHECK_INDEX(curr, ntokens);
    parName[1] = Tokenizer_token(tkz, curr++);

    if(curr < ntokens)
        eprintf("%s:%s:%d: extra token \"%s\" at end of line\n",
                 __FILE__,__func__,__LINE__, Tokenizer_token(tkz,curr));

    assert(strlen(childName) > 0);
    PopNode *childNode = PopNodeTab_get(poptbl, childName);
    if(childNode == NULL) {
        eprintf("%s:%s:%d: child segment \"%s\" undefined\n",
                 __FILE__,__func__,__LINE__, childName);
    }

    assert(strlen(parName[0]) > 0);
    PopNode *parNode0 = PopNodeTab_get(poptbl, parName[0]);
    if(parNode0 == NULL)
        eprintf("%s:%s:%d: parent segment \"%s\" undefined\n",
                 __FILE__,__func__,__LINE__, parName[0]);

    assert(strlen(parName[1]) > 0);
    PopNode *parNode1 = PopNodeTab_get(poptbl, parName[1]);
    if(parNode1 == NULL)
        eprintf("%s:%s:%d: parent segment \"%s\" undefined\n",
                 __FILE__,__func__,__LINE__, parName[1]);

    PopNode_mix(childNode, mPtr, mfree, parNode1, parNode0);
}

PopNode    *mktree(FILE * fp, SampNdx *sndx, LblNdx *lndx, ParStore *parstore,
                   ExoPar *exopar, Bounds *bnd, NodeStore *ns) {
    int         ntokens;
    char        buff[500];
    Tokenizer  *tkz = Tokenizer_new(50);

    PopNodeTab *poptbl = PopNodeTab_new();

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

        // Tokenize. Fields are separated by spaces, tabs, or "=".
        Tokenizer_split(tkz, buff, " \t=");

        ntokens = Tokenizer_strip(tkz, " \t=\n");
        if(ntokens == 0)
            continue;

        char *tok = Tokenizer_token(tkz, 0);
		if(0 == strcmp(tok, "twoN"))
			parseParam(tkz, TwoN, parstore, bnd, exopar, false);
		else if(0 == strcmp(tok, "time"))
			parseParam(tkz, Time, parstore, bnd, exopar, true);
		else if(0 == strcmp(tok, "mixFrac"))
			parseParam(tkz, MixFrac, parstore, bnd, exopar, false);
        else if(0 == strcmp(tok, "segment"))
            parseSegment(tkz, poptbl, sndx, lndx, parstore, ns);
        else if(0 == strcmp(tok, "mix"))
            parseMix(tkz, poptbl, parstore);
        else if(0 == strcmp(tok, "derive"))
            parseDerive(tkz, poptbl);
        else
            ILLEGAL_INPUT(tok);
    }

    // No more additions allowed to exopar.
    ExoPar_freeze(exopar);

    // Make sure the tree of populations has a single root. This
    // code iterates through all the nodes in the PopNodeTab, and
    // searches from each node back to the root. If all is well,
    // these searches all find the same root. Otherwise, it aborts
    // with an error.
    PopNodeTab_sanityCheck(poptbl, __FILE__, __LINE__);
    PopNode *root = PopNodeTab_root(poptbl);
    Tokenizer_free(tkz);
    PopNodeTab_free(poptbl);
    return root;
}

/// Count the number of "segment" statements in input file.
int countSegments(FILE * fp) {
    int         ntokens, nseg=0;
    char        buff[500];
    Tokenizer  *tkz = Tokenizer_new(50);

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
		if(0 == strcmp(tok, "segment"))
            ++nseg;
    }
    Tokenizer_free(tkz);
    return nseg;
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

    SampNdx sndx;
    SampNdx_init(&sndx);
	LblNdx lndx;
	LblNdx_init(&lndx);
	ParStore *parstore = ParStore_new();  // parameters
	Bounds   bnd = {
		.lo_twoN = 0.0,
		.hi_twoN = 1e7,
		.lo_t = 0.0,
		.hi_t = HUGE_VAL
	};

    int nseg = countSegments(fp);
    assert(6 == nseg);
    unitTstResult("countSegments", "OK");

    rewind(fp);

    PopNode  nodeVec[nseg];
    PopNode *root=NULL;
    ExoPar *ep = ExoPar_new();

    {
        NodeStore *ns = NodeStore_new(nseg, nodeVec);
        root = mktree(fp, &sndx, &lndx, parstore, ep, &bnd, ns);
        assert(root != NULL);
        NodeStore_free(ns);
    }

    if(verbose) {
        PopNode_print(stdout, root, 0);
        unsigned i;
        for(i=0; i < LblNdx_size(&lndx); ++i)
            printf("%2u %s\n", i, LblNdx_lbl(&lndx, i));
		printf("Used %d fixed parameters in \"parstore\".\n",
               ParStore_nFixed(parstore));
		printf("Used %d free parameters in \"parstore\".\n",
               ParStore_nFree(parstore));
    }

    unitTstResult("mktree", "needs more testing");

	ParStore_free(parstore);
    fclose(fp);
    unlink(tstFname);
    ExoPar_free(ep);
    return 0;
}
#endif
