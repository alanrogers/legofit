/**
 * @file parse.c
 * @brief Parse a .lgo file
 * @author Alan R. Rogers
 *
 * Consider the following tree of populations:
 *
 *      a-------|
 *              |ab--|
 *      b--|bb--|    |
 *         |         |abc--
 *         |c--------|
 *
 *  t = 0  1    3    5.5     inf
 *
 * In this tree, a, b, c, bb, ab, and abc represent segments of the
 * population tree.  The input file begins with a series of "segment"
 * statements that define each of the segments in the tree. The
 * segment statements also provide the time (backwards from the
 * present) at which the segment starts, and the size, twoN, of the
 * population. Optionally, it also provides the number of haploid
 * samples observed in this segment of the tree.
 *
 * The statements that follow the segment statements describe how the
 * segments are connected. The "mix" statement is used when a segment
 * originates as a mixture of two ancestral segments. The "derive"
 * statement is used when a segment derives from a single ancestral
 * segment.
 *
 * No segment can have more than two "parents" or more than two
 * "children".
 *
 * Here is input that would generate the tree above:
 *
 * time fixed T0=0
 * time free  Tc=1
 * time free  Tab=3
 * time free  Tabc=5.5
 * twoN free  twoNa=100
 * twoN fixed twoNb=123
 * twoN free  twoNc=213.4
 * twoN fixed twoNbb=32.1
 * twoN free  twoNab=222
 * twoN fixed twoNabc=1.2e2
 * mixFrac free Mc=0.8
 * segment a   t=T0     twoN=twoNa    samples=1
 * segment b   t=T0     twoN=twoNb    samples=2
 * segment c   t=Tc     twoN=twoNc    samples=1
 * segment bb  t=Tc     twoN=twoNbb
 * segment ab  t=Tab    twoN=twoNab
 * segment abc t=Tabc   twoN=twoNabc
 * mix    b  from bb + Mc * c
 * derive a  from ab
 * derive bb from ab
 * derive ab from abc
 * derive c  from abc
 *
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "popnodetab.h"
#include "lblndx.h"
#include "misc.h"
#include "parse.h"
#include "parstore.h"
#include "popnode.h"
#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/// Abort if token is missing
# define CHECK_TOKEN(tok, orig) {                                  \
        if((tok) == NULL) {                                        \
            fprintf(stderr, "%s:%d:"                               \
                    " input line incomplete in .lgo file.\n",      \
                    __FILE__,__LINE__);                            \
            fprintf(stderr,"  input: %s\n", (orig));               \
            exit(EXIT_FAILURE);                                    \
        }                                                          \
    }while(0);

/// Abort with an error message about illegal input.
#define ILLEGAL_INPUT(x,orig) do{                               \
        fprintf(stderr,"%s:%d: Illegal input: \"%s\"\n",        \
                __FILE__,__LINE__, (x));                        \
        fprintf(stderr,"  input: %s\n", (orig));                \
        exit(EXIT_FAILURE);                                     \
    }while(0)

int         getDbl(double *x, char **next, const char *orig);
int         getULong(unsigned long *x, char **next, const char *orig);
void		parseParam(char *next, enum ParamType type,
					   ParStore *parstore, Bounds *bnd,
                       const char *orig);
void        parseSegment(char *next, PopNodeTab *poptbl, SampNdx *sndx,
						 LblNdx *lndx, ParStore *parstore,
                         NodeStore *ns, const char *orig);
void        parseDerive(char *next, PopNodeTab *poptbl, const char *orig);
void        parseMix(char *next, PopNodeTab *poptbl, ParStore *parstore,
                     const char *orig);
int         get_one_line(size_t n, char buff[n], FILE *fp);

/// Interpret token i as a double.
/// @param[out] x points to variable into which double value will be
/// placed
/// @param[inout] next points to unparsed portion of input line
int getDbl(double *x, char **next, const char *orig) {
    char *tok, *end=NULL;
    tok = nextWhitesepToken(next);
    CHECK_TOKEN(tok, orig);
    *x = strtod(tok, &end);
    if(end!=tok && end != NULL && *end == '\0')
        return 0;               // success
    return 1;                   // failure
}

/// Interpret token i as an unsigned long integer.
/// @param[out] x points to variable into which value will be
/// placed
/// @param[inout] next points to unparsed portion of input line
/// integer.
int getULong(unsigned long *x, char **next, const char *orig) {
    char *tok, *end=NULL;
    tok = nextWhitesepToken(next);
    CHECK_TOKEN(tok, orig);
    *x = strtoul(tok, &end, 10);
    if(end != NULL && *end == '\0')
        return 0;               // success
    return 1;                   // failure
}

/// Parse a line of input defining a parameter
/// @param[inout] next points to unparsed portion of input line
/// @param[in] type the type of the current parameter: TwoN, Time, or MixFrac
/// @param[out] parstore structure that maintains info about
/// parameters
/// @param[in] bnd the bounds of each type of parameter
/// @parem[in] orig original input line
void		parseParam(char *next, enum ParamType type,
					   ParStore *parstore, Bounds *bnd,
                       const char *orig) {
    enum ParamStatus pstat = Free;

	// Read status of parameter
	{
		char *tok = nextWhitesepToken(&next);
        CHECK_TOKEN(tok, orig);
		if(0 == strcmp("fixed", tok))
			pstat = Fixed;
		else if(0 == strcmp("free", tok))
			pstat = Free;
        else if(0 == strcmp("gaussian", tok))
            pstat = Gaussian;
        else if(0 == strcmp("constrained", tok))
            pstat = Constrained;
		else {
			fprintf(stderr, "%s:%s:%d: got %s when expecting \"fixed\","
                    " \"free\", \"gaussian\", or \"constrained\".\n",
					__FILE__,__func__,__LINE__, tok);
            fprintf(stderr,"input: %s\n", orig);
			exit(EXIT_FAILURE);
		}
	}

	// Read parameter name, delimited by '='
	char *name = strsep(&next, "=");
    CHECK_TOKEN(name, orig);
    name = stripWhiteSpace(name);
    int i, ok=1;
	if( !isalpha(name[0]) )
        ok = 0;
    for(i=1; ok && name[i]!='\0'; ++i) {
        int c = name[i];
        if( !(isalnum(c) || c=='_') )
            ok = 0;
    }
    if(!ok) {
		fprintf(stderr,"%s:%d: \"%s\" is not a legal parameter name.\n",
				__FILE__,__LINE__, name);
        fprintf(stderr," Legal names consist of a letter followed by"
                " letters, digits, and underscores.\n");
        fprintf(stderr," Input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    char *formula;
    double value, sd = 0.0;
    if(pstat == Constrained) {
        formula = stripWhiteSpace(next);
        assert(formula != NULL);
    }else{
        // Read parameter value
        if(getDbl(&value, &next, orig)) {
            fflush(stdout);
            fprintf(stderr,
                    "%s:%s:%d:Can't parse double.\n",
                    __FILE__,__func__,__LINE__);
            fprintf(stderr,"input: %s\n", orig);
            exit(EXIT_FAILURE);
        }

        if(pstat == Gaussian) {

            // Read string "sd"
            char *sdstr = strsep(&next, "=");
            CHECK_TOKEN(sdstr, orig);
            sdstr = stripWhiteSpace(sdstr);
            if( 0 != strcmp(sdstr, "sd") ) {
                fprintf(stderr,"%s:%s:%d: \"%s\" should be \"sd\".\n",
                        __FILE__,__func__,__LINE__, sdstr);
                fprintf(stderr,"input: %s\n", orig);
                exit(EXIT_FAILURE);
            }

            // Read sd value
            if(getDbl(&sd, &next, orig)) {
                fflush(stdout);
                fprintf(stderr,
                        "%s:%s:%d:Can't parse double.\n",
                        __FILE__,__func__,__LINE__);
                fprintf(stderr,"input: %s\n", orig);
                exit(EXIT_FAILURE);
            }
            if(sd <= 0.0) {
                fprintf(stderr,"%s:%s:%d: Error sd=%lg <= 0.0\n",
                        __FILE__,__func__,__LINE__, sd);
                fprintf(stderr,"input: %s\n", orig);
                exit(EXIT_FAILURE);
            }
        }
    }

	// Allocate and initialize parameter in ParStore
    switch(pstat) {
    case Fixed:
		ParStore_addFixedPar(parstore, value, name);
        break;
    case Gaussian:
		ParStore_addGaussianPar(parstore, value, sd, name);
        break;
    case Constrained:
        ParStore_addConstrainedPar(parstore, formula, name);
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
                DIE("This shouldn't happen");
            }
            ParStore_addFreePar(parstore, value, lo, hi, name);
        }
        break;
    default:
        DIE("This shouldn't happen");
    }
}

/// Parse a line describing a segment of the population tree
/// @param[inout] next pointer to unparsed portion of input line
/// @param[inout] poptbl associates names of segments
/// with pointers to them.
/// @param[inout] sndx associates the index of each
/// sample with the PopNode object to which it belongs.
/// @param[inout] lndx associated index of each sample with its name
/// @param[out] parstore structure that maintains info about
/// parameters
/// @param[inout] ns allocates PopNode objects
void parseSegment(char *next, PopNodeTab *poptbl, SampNdx *sndx,
				  LblNdx *lndx, ParStore *parstore, NodeStore *ns,
                       const char *orig) {
    char *popName, *tok;
    double *tPtr, *twoNptr;
	ParamStatus tstat, twoNstat;
    unsigned long nsamples=0;

    // Read name of segment
    popName = nextWhitesepToken(&next);
    CHECK_TOKEN(popName, orig);

    // Read t
    tok = strsep(&next, "=");
    CHECK_TOKEN(tok, orig);
    tok = stripWhiteSpace(tok);
    if(0 != strcmp("t", tok)) {
        fprintf(stderr, "Got %s when expecting \"t\" on input:\n",
                tok);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }
    tok = nextWhitesepToken(&next);
    CHECK_TOKEN(tok, orig);
    tPtr = ParStore_findPtr(parstore, &tstat, tok);
	if(NULL == tPtr) {
		fprintf(stderr,"%s:%s:%d: Parameter \"%s\" is undefined\n",
				__FILE__,__func__,__LINE__,tok);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    // Read twoN
    tok = strsep(&next, "=");
    CHECK_TOKEN(tok, orig);
    tok = stripWhiteSpace(tok);
    if(0 != strcmp("twoN", tok)) {
        fprintf(stderr, "Got %s when expecting \"twoN\" on input:\n",
                tok);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }
    tok = nextWhitesepToken(&next);
    CHECK_TOKEN(tok, orig);
    tok = stripWhiteSpace(tok);
    twoNptr = ParStore_findPtr(parstore, &twoNstat, tok);
	if(NULL == twoNptr) {
		fprintf(stderr,"%s:%s:%dParameter \"%s\" is undefined\n",
				__FILE__,__func__,__LINE__, tok);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    // Read (optional) number of samples
    if(next) {
        tok = strsep(&next, "=");
        CHECK_TOKEN(tok, orig);
        tok = stripWhiteSpace(tok);
        if(0 != strcmp("samples", tok)) {
            fprintf(stderr, "%s:%s:%d: got %s when expecting \"samples\"\n",
                     __FILE__,__func__,__LINE__, tok);
            fprintf(stderr,"input: %s\n", orig);
            exit(EXIT_FAILURE);
        }
        if(getULong(&nsamples, &next, orig)) {
            fprintf(stderr, "%s:%s:%d: Can't parse unsigned int."
                    " Expecting value of \"samples\"\n",
                    __FILE__,__func__,__LINE__);
            fprintf(stderr,"input: %s\n", orig);
            exit(EXIT_FAILURE);
        }else {
            if(nsamples > MAXSAMP) {
                fprintf(stderr,
                        "%s:%s:%d: %lu samples is too many: max is %d:\n",
                         __FILE__,__func__,__LINE__, nsamples, MAXSAMP);
                fprintf(stderr,"input: %s\n", orig);
                exit(EXIT_FAILURE);
            }
        }
    }

    if(next){
        fprintf(stderr,"%s:%d: extra token(s) \"%s\" at end of line\n",
                 __FILE__,__LINE__, next);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    assert(strlen(popName) > 0);
    PopNode *thisNode = PopNode_new(twoNptr, twoNstat==Free,
                                    tPtr, tstat==Free, ns);
    if(0 != PopNodeTab_insert(poptbl, popName, thisNode)) {
        fprintf(stderr,"%s:%d: duplicate \"segment %s\"\n",
                 __FILE__,__LINE__, popName);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }
    LblNdx_addSamples(lndx, nsamples, popName);
    SampNdx_addSamples(sndx, nsamples, thisNode);
}

/// Parse a line of input describing a parent-offspring relationship
/// between two nodes.
/// @param[in] next unparsed portion of input line
/// @param[inout] poptbl associates names of segments
/// with pointers to them.
void parseDerive(char *next, PopNodeTab *poptbl,
                       const char *orig) {
    char *childName, *parName, *tok;

    // Read name of child
    childName = nextWhitesepToken(&next);
    CHECK_TOKEN(childName, orig);
    childName = stripWhiteSpace(childName);

    // Read "from"
    tok = nextWhitesepToken(&next);
    CHECK_TOKEN(tok, orig);
    if(0 != strcmp("from", tok)) {
        fprintf(stderr, "%s:%d: Got %s when expecting \"from\" on input:\n",
                __FILE__,__LINE__,tok);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    // Read name of parent
    parName = nextWhitesepToken(&next);
    CHECK_TOKEN(parName, orig);
    parName = stripWhiteSpace(parName);

    if(next) {
        fprintf(stderr,"%s:%d: extra tokens \"%s\" at end of line\n",
                 __FILE__,__LINE__, next);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    assert(strlen(childName) > 0);
    PopNode *childNode = PopNodeTab_get(poptbl, childName);
    if(childNode == NULL) {
        fprintf(stderr,"%s:%d: child segment \"%s\" undefined\n",
                 __FILE__,__LINE__, childName);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    assert(strlen(parName) > 0);
    PopNode *parNode = PopNodeTab_get(poptbl, parName);
    if(parNode == NULL) {
        fprintf(stderr, "%s:%d: parent segment \"%s\" undefined\n",
                 __FILE__,__LINE__, parName);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }
    PopNode_addChild(parNode, childNode);
}

/// Parse a line of input describing gene flow.
/// @param[inout] next unparsed portion of input line
/// @param[inout] poptbl associates names of segments
/// with pointers to them.
/// @param[out] parstore structure that maintains info about
/// parameters
void parseMix(char *next, PopNodeTab *poptbl, ParStore *parstore,
                       const char *orig) {
    char *childName, *parName[2], *tok;
    double *mPtr;
	ParamStatus mstat;

    // Read name of child
    childName = nextWhitesepToken(&next);
    CHECK_TOKEN(childName, orig);

    // Read word "from"
    tok = nextWhitesepToken(&next);
    CHECK_TOKEN(tok, orig);
    if(0 != strcmp("from", tok)) {
        fprintf(stderr,"%s:%d: got %s when expecting \"from\" on input:\n",
                 __FILE__,__LINE__, tok);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    // Read name of parent0
    parName[0] = strsep(&next, "+");
    CHECK_TOKEN(parName[0], orig);
    parName[0] = stripWhiteSpace(parName[0]);

    // Read mixture fraction, mPtr
    tok = strsep(&next, "*");
    CHECK_TOKEN(tok, orig);
    tok = stripWhiteSpace(tok);
    mPtr = ParStore_findPtr(parstore, &mstat, tok);
	if(NULL == mPtr) {
		fprintf(stderr,"%s:%s:%d: Parameter \"%s\" is undefined\n",
				__FILE__,__func__,__LINE__, tok);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    // Read name of parent1
    parName[1] = nextWhitesepToken(&next);
    CHECK_TOKEN(parName[1], orig);

    if(next) {
        fprintf(stderr, "%s:%d: extra token \"%s\" at end of line\n",
                 __FILE__,__LINE__, tok);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    assert(strlen(childName) > 0);
    PopNode *childNode = PopNodeTab_get(poptbl, childName);
    if(childNode == NULL) {
        fprintf(stderr, "%s:%d: child segment \"%s\" undefined\n",
                 __FILE__,__LINE__, childName);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    assert(strlen(parName[0]) > 0);
    PopNode *parNode0 = PopNodeTab_get(poptbl, parName[0]);
    if(parNode0 == NULL) {
        fprintf(stderr,"%s:%d: parent segment \"%s\" undefined\n",
                 __FILE__,__LINE__, parName[0]);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    assert(strlen(parName[1]) > 0);
    PopNode *parNode1 = PopNodeTab_get(poptbl, parName[1]);
    if(parNode1 == NULL) {
        fprintf(stderr,"%s:%d: parent segment \"%s\" undefined\n",
                 __FILE__,__LINE__, parName[1]);
        fprintf(stderr,"input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    PopNode_mix(childNode, mPtr, mstat==Free, parNode1, parNode0);
}

// Read a line into buff; strip comments and trailing whitespace.
// Return 0 on success; 1 on EOF.
int get_one_line(size_t n, char buff[n], FILE *fp) {
    if(fgets(buff, n, fp) == NULL)
        return 1;

    if(!strchr(buff, '\n') && !feof(fp)) {
        fprintf(stderr, "%s:%d: buffer overflow. buff size: %zu\n",
                __FILE__, __LINE__, n);
        fprintf(stderr,"input: %s\n", buff);
        exit(EXIT_FAILURE);
    }

    // strip trailing comments
    char *s = strchr(buff, '#');
    if(s)
        *s = '\0';

    // strip trailing whitespace
    for(s=buff; *s != '\0'; ++s)
        ;
    while(s > buff && isspace( *(s-1) ))
        --s;
    *s = '\0';

    return 0;
}

/// Parse an input file in .lgo format
/// @param[inout] fp input file pointer
/// @param[inout] sndx associates the index of each
/// sample with the PopNode object to which it belongs.
/// @param[inout] lndx associated index of each sample with its name
/// @param[out] parstore structure that maintains info about
/// parameters
/// @param[in] bnd the bounds of each type of parameter
/// @param[inout] ns allocates PopNode objects
PopNode    *mktree(FILE * fp, SampNdx *sndx, LblNdx *lndx, ParStore *parstore,
                   Bounds *bnd, NodeStore *ns) {
    char        orig[500], buff[500], buff2[500];
    char        *token, *next;

    PopNodeTab *poptbl = PopNodeTab_new();

    while(1) {
        if(1 == get_one_line(sizeof(buff), buff, fp))
            break;

        char *plus, *end;

        // If line ends with "+", then append next line
        while(1){
            end = buff + strlen(buff);
            assert(end < buff + sizeof(buff));
            plus = strrchr(buff, '+');
            if(plus==NULL || 1+plus != end)
                break;
            // line ends with plus: append next line
            if(1 == get_one_line(sizeof(buff2), buff2, fp)) {
                fprintf(stderr,"%s:%d: unexpected end of file\n",
                        __FILE__,__LINE__);
                exit(EXIT_FAILURE);
            }
            if(strlen(buff) + strlen(buff2) >= sizeof(buff)) {
                fprintf(stderr,"%s:%d: "
                        "buffer overflow on continuation line\n",
                        __FILE__,__LINE__);
                exit(EXIT_FAILURE);
            }
            strcat(buff, buff2);
        }

        snprintf(orig, sizeof orig, "%s", buff);

        // Get first whitespace-separated token
        next = stripWhiteSpace(buff);
        token = nextWhitesepToken(&next);
        if(token == NULL)
            continue;

		if(0 == strcmp(token, "twoN"))
			parseParam(next, TwoN, parstore, bnd, orig);
		else if(0 == strcmp(token, "time"))
			parseParam(next, Time, parstore, bnd, orig);
		else if(0 == strcmp(token, "mixFrac"))
			parseParam(next, MixFrac, parstore, bnd, orig);
        else if(0 == strcmp(token, "segment"))
            parseSegment(next, poptbl, sndx, lndx, parstore, ns, orig);
        else if(0 == strcmp(token, "mix"))
            parseMix(next, poptbl, parstore, orig);
        else if(0 == strcmp(token, "derive"))
            parseDerive(next, poptbl, orig);
        else
            ILLEGAL_INPUT(token, orig);
    }

    // Make sure the tree of populations has a single root. This
    // code iterates through all the nodes in the PopNodeTab, and
    // searches from each node back to the root. If all is well,
    // these searches all find the same root. Otherwise, it aborts
    // with an error.
    PopNode *root = PopNodeTab_check_and_root(poptbl, __FILE__, __LINE__);
    PopNodeTab_free(poptbl);
    return root;
}

/// Count the number of "segment" statements in input file.
int countSegments(FILE * fp) {
    int         nseg=0;
    char        orig[500], buff[500];
    char        *tok, *next;

    while(1) {
        if(fgets(buff, sizeof(buff), fp) == NULL)
            break;
        next = buff;

        if(!strchr(buff, '\n') && !feof(fp)) {
            fprintf(stderr, "ERR@%s:%d: input buffer overflow."
                    " buff size: %zu\n", __FILE__, __LINE__, sizeof(buff));
            fprintf(stderr,"input: %s\n", orig);
            exit(EXIT_FAILURE);
        }

        // strip trailing comments
        char *comment = strchr(buff, '#');
        if(comment)
            *comment = '\0';

        snprintf(orig, sizeof orig, "%s", buff);
        tok = nextWhitesepToken(&next);
        if(tok==NULL)
            continue;

		if(0 == strcmp(tok, "segment"))
            ++nseg;
    }
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
    "twoN free   twoNa=100\n"
    "twoN fixed  twoNb=123\n"
    "twoN free   twoNc=213.4\n"
    "twoN fixed  twoNbb=32.1\n"
    "twoN constrained twoNab=100 + -1.2*Tab\n"
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

    // test getDbl
    fprintf(stderr,"%s:%s:%d\n", __FILE__,__func__,__LINE__);
    char buff[100], *next;
    double x;
    strcpy(buff, " +1.23 ");
    next = buff;
    assert(0==getDbl(&x, &next, buff));
    assert(Dbl_near(x, 1.23));
    strcpy(buff, " -1.23e-4 ");

    next=buff;
    assert(0==getDbl(&x, &next, buff));
    assert(Dbl_near(x, -1.23e-4));
    fprintf(stderr,"%s:%s:%d\n", __FILE__,__func__,__LINE__);
    unitTstResult("getDbl", "OK");

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
    fprintf(stderr,"nseg=%d\n", nseg);
    assert(6 == nseg);
    unitTstResult("countSegments", "OK");

    rewind(fp);

    PopNode  nodeVec[nseg];
    PopNode *root=NULL;

    {
        NodeStore *ns = NodeStore_new(nseg, nodeVec);
        root = mktree(fp, &sndx, &lndx, parstore, &bnd, ns);
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
		printf("Used %d Gaussian parameters in \"parstore\".\n",
               ParStore_nGaussian(parstore));
		printf("Used %d constrained parameters in \"parstore\".\n",
               ParStore_nConstrained(parstore));
    }

    unitTstResult("mktree", "needs more testing");

	ParStore_free(parstore);
    fclose(fp);
    unlink(tstFname);
    return 0;
}
#endif
