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
 * present in generations) at which the segment starts, and the size,
 * twoN, of the  population. Optionally, it also provides the number
 * of haploid samples observed in this segment of the tree.
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

#include "error.h"
#include "lblndx.h"
#include "misc.h"
#include "network.h"
#include "param.h"
#include "parse.h"
#include "parstore.h"
#include "ptrqueue.h"
#include "sampndx.h"
#include "strptrmap.h"
#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Site pattern representing the union of all samples.
extern tipId_t union_all_samples;

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

/// Abort with an error message about duplicate parameter definition
#define DUPLICATE_PAR(x,orig) do{                                  \
        fprintf(stderr,"%s:%d: Duplicate parameter def: \"%s\"\n", \
                __FILE__,__LINE__, (x));                           \
        fprintf(stderr,"  input: %s\n", (orig));                   \
        exit(EXIT_FAILURE);                                        \
    }while(0)

/// Input statements out of order
#define ORDER_ERROR(word) do{                                           \
        fprintf(stderr,"%s:%d: Order error at \"%s\" in .lgo file.\n",  \
                __FILE__,__LINE__, (word));                             \
        fprintf(stderr,"Parameter definitions should come first\n"      \
                "   (twoN, time, mixFrac, and param),\n"                \
                "   then \"segment\" statements, and finally \"mix\"\n" \
                "   and \"derive\" statements.\n");                     \
        exit(EXIT_FAILURE);                                             \
    }while(0)

int getDbl(double *x, char **next, const char *orig);
int getULong(unsigned long *x, char **next, const char *orig);
int getRange(double x[2], char **next, const char *orig);
void parseParam(char *next, unsigned type, StrPtrMap *parmap,
                PtrQueue *fixedQ, PtrQueue *freeQ, PtrQueue *constrQ,
                Bounds * bnd, const char *orig);
void parseSegment(char *next, StrPtrMap * popmap, SampNdx * sndx,
                  LblNdx * lndx, ParStore * parstore,
                  const char *orig);
void parseDerive(char *next, StrPtrMap * popmap, ParStore * parstore,
                 const char *orig);
void parseMix(char *next, StrPtrMap * popmap, ParStore * parstore,
              const char *orig);
static int get_one_line(size_t n, char buff[n], FILE * fp);

/// Interpret token i as a double.
/// @param[out] x points to variable into which double value will be
/// placed
/// @param[inout] next points to unparsed portion of input line
int getDbl(double *x, char **next, const char *orig) {
    char *tok, *end = NULL;
    tok = nextWhitesepToken(next);
    CHECK_TOKEN(tok, orig);
    *x = strtod(tok, &end);
    if(end != tok && end != NULL && *end == '\0')
        return 0;               // success
    return 1;                   // failure
}

/// Interpret token i as an unsigned long integer.
/// @param[out] x points to variable into which value will be
/// placed
/// @param[inout] next points to unparsed portion of input line
/// integer.
int getULong(unsigned long *x, char **next, const char *orig) {
    char *tok, *end = NULL;
    tok = nextWhitesepToken(next);
    CHECK_TOKEN(tok, orig);
    *x = strtoul(tok, &end, 10);
    if(end != NULL && *end == '\0')
        return 0;               // success
    return 1;                   // failure
}

/// Read a range in form "[ 12, 34 ]". Return 0 on success
/// or 1 if range is not present. Abort if first character is "["
/// but the rest of the string is not interpretable as a range.
int getRange(double x[2], char **next, const char *orig) {
    while(isspace(**next))
        *next += 1;

    if(**next != '[')
        return 1;  // no range

    *next += 1;
    char *tok, *end=NULL;

    // Read lower bound
    tok = strsep(next, ",");
    CHECK_TOKEN(tok, orig);
    x[0] = strtod(tok, &end);
    if(end == tok || end == NULL || *end != '\0') {
        fprintf(stderr,"%s:%d: error reading lower bound.\n",
                __FILE__,__LINE__);
        fprintf(stderr,"  input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    // Read upper bound
    tok = strsep(next, "]");
    CHECK_TOKEN(tok, orig);
    tok = stripWhiteSpace(tok);
    x[1] = strtod(tok, &end);
    if(end == tok || end == NULL || *end != '\0') {
        fprintf(stderr,"%s:%d: error reading upper bound.\n",
                __FILE__,__LINE__);
        fprintf(stderr,"  input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    return 0;  // success
}

/// Parse a line of input defining a parameter
/// @param[inout] next points to unparsed portion of input line
/// @param[in] ptype TWON, TIME, or MIXFRAC
/// @param[out] parstore structure that maintains info about
/// parameters
/// @param[in] bnd the bounds of each type of parameter
/// @param[in] orig original input line
void parseParam(char *next, unsigned ptype, StrPtrMap *parmap,
                PtrQueue *fixedQ, PtrQueue *freeQ, PtrQueue *constrQ,
                Bounds * bnd, const char *orig) {
    // Read type of parameter
    {
        char *tok = nextWhitesepToken(&next);
        CHECK_TOKEN(tok, orig);
        if(0 == strcmp("fixed", tok))
            ptype |= FIXED;
        else if(0 == strcmp("free", tok))
            ptype |= FREE;
        else if(0 == strcmp("constrained", tok))
            ptype |= CONSTRAINED;
        else {
            fprintf(stderr, "%s:%s:%d: got %s when expecting \"fixed\","
                    " \"free\", or \"constrained\".\n",
                    __FILE__, __func__, __LINE__, tok);
            fprintf(stderr, "input: %s\n", orig);
            exit(EXIT_FAILURE);
        }
    }

    double range[2];
    int status = getRange(range, &next, orig);
    int gotRange = (status ? 0 : 1);

    if(gotRange && !(ptype & FREE)) {
        fprintf(stderr,"%s:%d: only free variables can have ranges\n",
                __FILE__,__LINE__);
        fprintf(stderr,"  input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    // Read parameter name, delimited by '='
    char *name = strsep(&next, "=");
    CHECK_TOKEN(name, orig);
    name = stripWhiteSpace(name);
    if( !legalName(name) ) {
        fprintf(stderr, "%s:%d: \"%s\" is not a legal parameter name.\n",
                __FILE__, __LINE__, name);
        fprintf(stderr, " Legal names consist of a letter followed by"
                " letters, digits, underscores, colons, or periods.\n");
        fprintf(stderr, " Input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    char *formula;
    double value=0.0;
    if(ptype & CONSTRAINED) {
        formula = stripWhiteSpace(next);
        assert(formula != NULL);
    } else {
        // Read parameter value
        if(getDbl(&value, &next, orig)) {
            fflush(stdout);
            fprintf(stderr,
                    "%s:%s:%d:Can't parse double.\n",
                    __FILE__, __func__, __LINE__);
            fprintf(stderr, "input: %s\n", orig);
            exit(EXIT_FAILURE);
        }
    }

    // Allocate and initialize parameter in queue
    Param *par;
    if(ptype & FIXED) {
        par = Param_new(name, value, value, value, ptype, NULL);
        PtrQueue_push(fixedQ, par);
        status = StrPtrMap_insert(parmap, name, par);
        if(status)
            DUPLICATE_PAR(name, orig);
    }else if(ptype & CONSTRAINED) {
        par = Param_new(name, value, -DBL_MAX, DBL_MAX, ptype, formula);
        PtrQueue_push(constrQ, par);
        status = StrPtrMap_insert(parmap, name, par);
        if(status)
            DUPLICATE_PAR(name, orig);
    }else if(ptype & FREE) {
        double lo, hi;
        if(gotRange) {
            lo = range[0];
            hi = range[1];
        }else{
            if(ptype & TWON) {
                lo = bnd->lo_twoN;
                hi = bnd->hi_twoN;
            }else if(ptype & TIME) {
                lo = bnd->lo_t;
                hi = bnd->hi_t;
            }else if (ptype & MIXFRAC) {
                lo = 0.0;
                hi = 1.0;
            }else if (ptype & ARBITRARY) {
                lo = -DBL_MAX;
                hi = DBL_MAX;
            }else
                DIE("This shouldn't happen");
        }
        par = Param_new(name, value, lo, hi, ptype, NULL);
        PtrQueue_push(freeQ, par);
        status = StrPtrMap_insert(parmap, name, par);
        if(status)
            DUPLICATE_PAR(name, orig);
    }else
        DIE("This shouldn't happen");
}

/// Parse a line describing a segment of the population tree
/// @param[inout] next pointer to unparsed portion of input line
/// @param[inout] popmap associates names of segments
/// with pointers to them.
/// @param[inout] sndx associates the index of each
/// sample with the node to which it belongs.
/// @param[inout] lndx associated index of each sample with its name
/// @param[out] parstore structure that maintains info about
/// parameters
void parseSegment(char *next, StrPtrMap * popmap, SampNdx * sndx,
                  LblNdx * lndx, ParStore * parstore, const char *orig) {
    char *popName, *tok;
    int start_i, twoN_i;
    unsigned long nsamples = 0;

    // Read name of segment
    popName = nextWhitesepToken(&next);
    CHECK_TOKEN(popName, orig);

    // Read t
    tok = strsep(&next, "=");
    CHECK_TOKEN(tok, orig);
    tok = stripWhiteSpace(tok);
    if(0 != strcmp("t", tok)) {
        fprintf(stderr, "Got %s when expecting \"t\" on input:\n", tok);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }
    tok = nextWhitesepToken(&next);
    CHECK_TOKEN(tok, orig);
    start_i = ParStore_getIndex(parstore, tok);
    if(start_i < 0) {
        fprintf(stderr, "%s:%s:%d: Parameter \"%s\" is undefined\n",
                __FILE__, __func__, __LINE__, tok);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }
    // Read twoN
    tok = strsep(&next, "=");
    CHECK_TOKEN(tok, orig);
    tok = stripWhiteSpace(tok);
    if(0 != strcmp("twoN", tok)) {
        fprintf(stderr, "Got %s when expecting \"twoN\" on input:\n", tok);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }
    tok = nextWhitesepToken(&next);
    CHECK_TOKEN(tok, orig);
    tok = stripWhiteSpace(tok);
    twoN_i = ParStore_getIndex(parstore, tok);
    if(twoN_i < 0) {
        fprintf(stderr, "%s:%s:%dParameter \"%s\" is undefined\n",
                __FILE__, __func__, __LINE__, tok);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    // Read (optional) number of samples
    if(next) {
        tok = strsep(&next, "=");
        CHECK_TOKEN(tok, orig);
        tok = stripWhiteSpace(tok);
        if(0 != strcmp("samples", tok)) {
            fprintf(stderr, "%s:%s:%d: got %s when expecting \"samples\"\n",
                    __FILE__, __func__, __LINE__, tok);
            fprintf(stderr, "input: %s\n", orig);
            exit(EXIT_FAILURE);
        }
        if(getULong(&nsamples, &next, orig)) {
            fprintf(stderr, "%s:%s:%d: Can't parse unsigned int."
                    " Expecting value of \"samples\"\n",
                    __FILE__, __func__, __LINE__);
            fprintf(stderr, "input: %s\n", orig);
            exit(EXIT_FAILURE);
        } else {
            if(nsamples > MAXSAMP) {
                fprintf(stderr,
                        "%s:%s:%d: %lu samples is too many: max is %d:\n",
                        __FILE__, __func__, __LINE__, nsamples, MAXSAMP);
                fprintf(stderr, "input: %s\n", orig);
                exit(EXIT_FAILURE);
            }
        }
    }

    if(next) {
        fprintf(stderr, "%s:%d: extra token(s) \"%s\" at end of line\n",
                __FILE__, __LINE__, next);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    assert(strlen(popName) > 0);
    void *thisNode = Node_new(twoN_i, start_i, parstore, popName);
    if(0 != StrPtrMap_insert(popmap, popName, thisNode)) {
        fprintf(stderr, "%s:%d: duplicate \"segment %s\"\n",
                __FILE__, __LINE__, popName);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }
    LblNdx_addSamples(lndx, nsamples, popName);
    SampNdx_addSamples(sndx, nsamples, thisNode);
}

/// Parse a line of input describing a parent-offspring relationship
/// between two nodes.
/// @param[in] next unparsed portion of input line
/// @param[inout] popmap associates names of segments
/// with pointers to them.
void parseDerive(char *next, StrPtrMap * popmap, ParStore * parstore,
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
                __FILE__, __LINE__, tok);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }
    // Read name of parent
    parName = nextWhitesepToken(&next);
    CHECK_TOKEN(parName, orig);
    parName = stripWhiteSpace(parName);

    if(next) {
        fprintf(stderr, "%s:%d: extra tokens \"%s\" at end of line\n",
                __FILE__, __LINE__, next);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    assert(strlen(childName) > 0);
    void *childNode = StrPtrMap_get(popmap, childName);
    if(childNode == NULL) {
        fprintf(stderr, "%s:%d: child segment \"%s\" undefined\n",
                __FILE__, __LINE__, childName);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    assert(strlen(parName) > 0);
    void *parNode = StrPtrMap_get(popmap, parName);
    if(parNode == NULL) {
        fprintf(stderr, "%s:%d: parent segment \"%s\" undefined\n",
                __FILE__, __LINE__, parName);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }
    int status = Node_addChild(parNode, childNode);
    if(status) {
        char buff[200];
        mystrerror_r(status, buff, sizeof(buff));
        fprintf(stderr,"%s:%d: %s\n", __FILE__,__LINE__,buff);
        exit(EXIT_FAILURE);
    }
}

/// Parse a line of input describing gene flow.
/// @param[inout] next unparsed portion of input line
/// @param[inout] popmap associates names of segments
/// with pointers to them.
/// @param[out] parstore structure that maintains info about
/// parameters
void parseMix(char *next, StrPtrMap * popmap, ParStore * parstore,
              const char *orig) {
    char *childName, *parName[2], *tok;
    int mix_i;

    // Read name of child
    childName = nextWhitesepToken(&next);
    CHECK_TOKEN(childName, orig);

    // Read word "from"
    tok = nextWhitesepToken(&next);
    CHECK_TOKEN(tok, orig);
    if(0 != strcmp("from", tok)) {
        fprintf(stderr, "%s:%d: got %s when expecting \"from\" on input:\n",
                __FILE__, __LINE__, tok);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }
    // Read name of parent0
    parName[0] = strsep(&next, "+");
    CHECK_TOKEN(parName[0], orig);
    parName[0] = stripWhiteSpace(parName[0]);

    // Read mixture fraction
    tok = strsep(&next, "*");
    CHECK_TOKEN(tok, orig);
    tok = stripWhiteSpace(tok);
    mix_i = ParStore_getIndex(parstore, tok);
    if(mix_i < 0) {
        fprintf(stderr, "%s:%s:%d: Parameter \"%s\" is undefined\n",
                __FILE__, __func__, __LINE__, tok);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }
    // Read name of parent1
    parName[1] = nextWhitesepToken(&next);
    CHECK_TOKEN(parName[1], orig);

    if(next) {
        fprintf(stderr, "%s:%d: extra token \"%s\" at end of line\n",
                __FILE__, __LINE__, tok);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    assert(strlen(childName) > 0);
    void *childNode = StrPtrMap_get(popmap, childName);
    if(childNode == NULL) {
        fprintf(stderr, "%s:%d: child segment \"%s\" undefined\n",
                __FILE__, __LINE__, childName);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    assert(strlen(parName[0]) > 0);
    void *parNode0 = StrPtrMap_get(popmap, parName[0]);
    if(parNode0 == NULL) {
        fprintf(stderr, "%s:%d: parent segment \"%s\" undefined\n",
                __FILE__, __LINE__, parName[0]);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    assert(strlen(parName[1]) > 0);
    void *parNode1 = StrPtrMap_get(popmap, parName[1]);
    if(parNode1 == NULL) {
        fprintf(stderr, "%s:%d: parent segment \"%s\" undefined\n",
                __FILE__, __LINE__, parName[1]);
        fprintf(stderr, "input: %s\n", orig);
        exit(EXIT_FAILURE);
    }

    int status = Node_mix(childNode, mix_i, parNode1, parNode0,
                          parstore);
    switch(status){
    case 0:
        break;
    case TOO_MANY_PARENTS:
        fprintf(stderr,"%s:%d: node has too many parents\n",
                __FILE__,__LINE__);
        fprintf(stderr, " Input: %s\n", orig);
        exit(EXIT_FAILURE);
    case TOO_MANY_CHILDREN:
        fprintf(stderr,"%s:%d: node has too many children\n",
                __FILE__,__LINE__);
        fprintf(stderr, " Input: %s\n", orig);
        exit(EXIT_FAILURE);
    case DATE_MISMATCH:
        fprintf(stderr,"%s:%d: date mismatch\n",__FILE__,__LINE__);
        fprintf(stderr, " Input: %s\n", orig);
        exit(EXIT_FAILURE);
    default:
        fprintf(stderr,"%s:%d: unknown error (%d)\n",
                __FILE__,__LINE__,status);
        fprintf(stderr, " Input: %s\n", orig);
        exit(EXIT_FAILURE);
    }
}

// Read a line into buff, skipping blank lines; strip comments and
// trailing whitespace.  Return 0 on success; EOF on end of file.
static int get_one_line(size_t n, char buff[n], FILE * fp) {
    do{
        if(fgets(buff, n, fp) == NULL)
            return EOF;

        if(!strchr(buff, '\n') && !feof(fp)) {
            fprintf(stderr, "%s:%d: buffer overflow. buff size: %zu\n",
                    __FILE__, __LINE__, n);
            fprintf(stderr, "input: %s\n", buff);
            exit(EXIT_FAILURE);
        }
        // strip trailing comments
        char *s = strchr(buff, '#');
        if(s)
            *s = '\0';

        // strip trailing whitespace
        for(s = buff; *s != '\0'; ++s) ;
        while(s > buff && isspace(*(s - 1)))
            --s;
        *s = '\0';
    }while(*buff == '\0');

    return 0;
}

/// Parse an input file in .lgo format
/// @param[inout] fp input file pointer
/// @param[inout] sndx associates the index of each
/// sample with the node to which it belongs.
/// @param[inout] lndx associated index of each sample with its name
/// @param[out] parstore structure that maintains info about
/// parameters
/// @param[in] bnd the bounds of each type of parameter
PtrPair mktree(FILE * fp, SampNdx * sndx, LblNdx * lndx, Bounds * bnd) {  
    char orig[2048], buff[2048], buff2[2048];
    char *token, *next;
    ParStore *parstore = NULL;

    StrPtrMap *popmap = StrPtrMap_new();
    StrPtrMap *parmap = StrPtrMap_new();

    // Queues for fixed parameters, free ones, and constrained ones.
    // Used during parsing. NULL after that.
    PtrQueue *fixedQ = PtrQueue_new();
    PtrQueue *freeQ = PtrQueue_new();
    PtrQueue *constrQ = PtrQueue_new();

    while(1) {
        if(EOF == get_one_line(sizeof(buff), buff, fp))
            break;

        // If line ends with a binary operator ("+-*/"), then append
        // next line.
        while(1) {
            char *end = buff + strlen(buff);
            // Check if last character is a binary operator.
            // No need to strip trailing whitespace, because that's
            // done in get_one_line.
            assert(end==buff || !isspace(*(end-1)));
            if(end>buff && strchr("+-*/", *(end-1)) == NULL)
                break;
            // line empty or ends with binary operator: append next line
            if(EOF == get_one_line(sizeof(buff2), buff2, fp)) {
                fprintf(stderr, "%s:%d: unexpected end of file\n",
                        __FILE__, __LINE__);
                fprintf(stderr,"  prev line: \"%s\"\n", buff);
                exit(EXIT_FAILURE);
            }
            if(strlen(buff) + 1 + strlen(buff2) >= sizeof(buff)) {
                fprintf(stderr, "%s:%d: "
                        "buffer overflow on continuation line\n",
                        __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            strcat(buff, " "); // add space after operator
            strcat(buff, buff2);
        }

        collapse_whitespace(buff);

        snprintf(orig, sizeof orig, "%s", buff);

        // Get first whitespace-separated token
        next = stripWhiteSpace(buff);
        token = nextWhitesepToken(&next);
        if(token == NULL)
            continue;

        if(0 == strcmp(token, "twoN")) {
            if(parstore != NULL)
               ORDER_ERROR("twoN");
            parseParam(next, TWON, parmap, fixedQ, freeQ, constrQ, bnd, orig);
        }else if(0 == strcmp(token, "time")) {
            if(parstore != NULL)
               ORDER_ERROR("time");
            parseParam(next, TIME, parmap, fixedQ, freeQ, constrQ, bnd, orig);
        }else if(0 == strcmp(token, "mixFrac")) {
            if(parstore != NULL)
               ORDER_ERROR("mixFrac");
            parseParam(next, MIXFRAC, parmap, fixedQ, freeQ, constrQ, bnd,
                       orig); 
        }else if(0 == strcmp(token, "param")) {
            if(parstore != NULL)
               ORDER_ERROR("param");
            parseParam(next, ARBITRARY, parmap, fixedQ, freeQ, constrQ, bnd,
                       orig);
        }else if(0 == strcmp(token, "segment")) {
            if(parstore==NULL) {
                // 1st segment: allocate parstore and free
                // fixedQ, freeQ, and constrQ.
                assert(fixedQ != NULL);
                parstore = ParStore_new(fixedQ, freeQ, constrQ);
                CHECKMEM(parstore);
                PtrQueue_free(fixedQ);
                PtrQueue_free(freeQ);
                PtrQueue_free(constrQ);
                fixedQ = freeQ = constrQ = NULL;
            }
            parseSegment(next, popmap, sndx, lndx, parstore, orig);
        } else if(0 == strcmp(token, "mix")) {
            if(parstore == NULL)
               ORDER_ERROR("mix");
            parseMix(next, popmap, parstore, orig);
        }else if(0 == strcmp(token, "derive")) {
            if(parstore == NULL)
               ORDER_ERROR("derive");
            parseDerive(next, popmap, parstore, orig);
        }
        else
            ILLEGAL_INPUT(token, orig);
    }

    // Make sure the tree of populations has a single root. This code
    // iterates through all the nodes, and searches from each node
    // back to the root. If all is well, these searches all find the
    // same root. Otherwise, abort with an error.
    void *root = NULL;
    {
        /// Check the sanity of each node and make sure there is only one
        /// root.
        unsigned long n = StrPtrMap_size(popmap);
        void *nodes[n];
        void *curr;
        StrPtrMap_ptrArray(popmap, n, nodes);
        for(unsigned long j=0; j < n; ++j) {
            curr = Node_root(nodes[j]);
            if(root == NULL)
                root = curr;
            else if (root != curr) {
                fprintf(stderr,
                        "%s:%d: Pop tree has multiple roots.\n",
                        __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        }
    }

    // site pattern representing the union of all samples
    unsigned nsamples = SampNdx_size(sndx);
    union_all_samples = (1u << nsamples) - 1;
    
    StrPtrMap_free(popmap);
    StrPtrMap_free(parmap);

    assert(root);
    assert(parstore);
    PtrPair ptrpair = {root, parstore};
    return ptrpair;
}

/// Count the number of "segment" statements in input file.
int countSegments(FILE * fp) {
    int nseg = 0;
    char orig[500], buff[500];
    char *tok, *next;

    while(1) {
        if(fgets(buff, sizeof(buff), fp) == NULL)
            break;
        next = buff;

        if(!strchr(buff, '\n') && !feof(fp)) {
            fprintf(stderr, "ERR@%s:%d: input buffer overflow."
                    " buff size: %zu\n", __FILE__, __LINE__, sizeof(buff));
            fprintf(stderr, "input: %s\n", orig);
            exit(EXIT_FAILURE);
        }
        // strip trailing comments
        char *comment = strchr(buff, '#');
        if(comment)
            *comment = '\0';

        snprintf(orig, sizeof orig, "%s", buff);
        tok = nextWhitesepToken(&next);
        if(tok == NULL)
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

    int verbose = 0;

    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xparse [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    Network_init(STOCHASTIC);
    const char *tstFname = "mktree-tmp.lgo";
    FILE *fp = fopen(tstFname, "w");
    fputs(tstInput, fp);
    fclose(fp);
    fp = fopen(tstFname, "r");
    if(fp == NULL) {
        fprintf(stderr, "%s:%d: Can't open file \"%s\"\n",
                __FILE__, __LINE__, tstFname);
        exit(1);
    }

    // test getDbl
    char buff[100], *next;
    double x;
    strcpy(buff, " +1.23 ");
    next = buff;
    assert(0 == getDbl(&x, &next, buff));
    assert(Dbl_near(x, 1.23));
    strcpy(buff, " -1.23e-4 ");

    next = buff;
    assert(0 == getDbl(&x, &next, buff));
    assert(Dbl_near(x, -1.23e-4));
    unitTstResult("getDbl", "OK");

    // test getULong
    long unsigned lu;
    strcpy(buff, " +123 ");
    next = buff;
    assert(0 == getULong(&lu, &next, buff));
    assert(lu == 123);
    unitTstResult("getULong", "OK");

    // test getRange
    double y[2];
    strcpy(buff, " [ 12.34, 2.3e2 ] ");
    next = buff;
    assert(0 == getRange(y, &next, buff));
    assert(y[0] == 12.34);
    assert(y[1] == 2.3e2);

    strcpy(buff, "  12.34, 2.3e2 ] ");
    next = buff;
    assert(1 == getRange(y, &next, buff));
    assert(*next == '1');

    unitTstResult("getRange", "OK");

    SampNdx sndx;
    SampNdx_init(&sndx);
    LblNdx lndx;
    LblNdx_init(&lndx);
    Bounds bnd = {
        .lo_twoN = 0.0,
        .hi_twoN = 1e7,
        .lo_t = 0.0,
        .hi_t = HUGE_VAL
    };

    int nseg = countSegments(fp);
    assert(6 == nseg);
    unitTstResult("countSegments", "OK");

    rewind(fp);

    PopNode *root;
    ParStore *parstore;

    {
        PtrPair pp = mktree(fp, &sndx, &lndx, &bnd);
        root = pp.a;
        parstore = pp.b;
        assert(root != NULL);
        assert(parstore != NULL);
    }

    if(verbose) {
        Node_print(root, stdout, 0);
        unsigned i;
        for(i = 0; i < LblNdx_size(&lndx); ++i)
            printf("%2u %s\n", i, LblNdx_lbl(&lndx, i));
        printf("Used %d fixed parameters in \"parstore\".\n",
               ParStore_nFixed(parstore));
        printf("Used %d free parameters in \"parstore\".\n",
               ParStore_nFree(parstore));
        printf("Used %d constrained parameters in \"parstore\".\n",
               ParStore_nConstrained(parstore));
    }

    PopNode_free(root);
    assert(parstore);
    ParStore_free(parstore);

    unitTstResult("mktree", "OK");

    fclose(fp);
    unlink(tstFname);
    return 0;
}
#endif
