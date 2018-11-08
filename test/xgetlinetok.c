/**
 * @file xgetlinetok.c
 * @author Alan R. Rogers
 * @brief Test getlinetok.c.
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "getlinetok.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

const char *tstInput =
    "#chr  pos  ref  alt  raf\n"
    "1  1  a  t  0\n"
    "onefield\n";

int main(int argc, char **argv) {
    size_t      buffsize = 2;
    int         maxtokens = 0;
    int         verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0) 
            DIE("usage: xgetlinetok [-v]");
        verbose = 1;
        break;
    default:
        DIE("usage: xgetlinetok [-v]");
    }

    const char *fname = "getlinetok.tmp";
    FILE *fp = fopen(fname,"w");
    fputs(tstInput, fp);
    fclose(fp);
    fp = fopen(fname,"r");

    GetLineTok *glt = GetLineTok_new(buffsize, maxtokens, fp);
    assert(glt);

    assert(0 == GetLineTok_ntokens(glt));
    const char *sep = " ";
    const char *extraneous = "\n";
    
    if(verbose) {
        fprintf(stderr,"sep=\"%s\"\n", sep);
        fprintf(stderr,"extraneous=\"%s\"\n", extraneous);
    }

    assert(0 == GetLineTok_next(glt, sep, extraneous));
    assert(5 == GetLineTok_ntokens(glt));
    assert(0 == strcmp("raf", GetLineTok_token(glt, 4)));
    if(verbose)
        GetLineTok_print(glt, stderr);

    assert(0 == GetLineTok_next(glt, sep, extraneous));
    assert(5 == GetLineTok_ntokens(glt));
    assert(0 == strcmp("0", GetLineTok_token(glt, 4)));
    if(verbose)
        GetLineTok_print(glt, stderr);

    assert(0 == GetLineTok_next(glt, sep, extraneous));
    assert(1 == GetLineTok_ntokens(glt));
    assert(0 == strcmp("onefield", GetLineTok_token(glt, 0)));
    if(verbose)
        GetLineTok_print(glt, stderr);

    assert(EOF == GetLineTok_next(glt, sep, extraneous));
    assert(0 == GetLineTok_ntokens(glt));

    GetLineTok_free(glt);
    fclose(fp);
    remove(fname);
    unitTstResult("GetLineTok", "OK");

    return 0;
}
