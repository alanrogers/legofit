/**
 * @file xtokenizer.c
 * @author Alan R. Rogers
 * @brief Test tokenizer.c.
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "tokenizer.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {
    int         maxtokens = 1;
    int         ntokens, verbose = 0, i;
    Tokenizer  *tkz;
    const char *sep;
    char        str[100];

#ifdef NDEBUG
    eprintf("ERR@%s:%d:"
            "Unit tests must be compiled without -DNDEBUG flag\n",
            __FILE__, __LINE__);
#endif

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xtokenizer [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xtokenizer [-v]\n");
    }

    tkz = Tokenizer_new(maxtokens);
    CHECKMEM(tkz);

    int ntok = Tokenizer_ntokens(tkz);
    fprintf(stderr,"ntok=%d\n", ntok);

    assert(Tokenizer_ntokens(tkz) == 0);
    strcpy(str, "; now;   \t: \t is       ::the  ,time  \n,");
    sep = ";:,";
    if(verbose) {
        printf("sep=\"%s\"\n", sep);
        printf("str=\"%s\"\n", str);
    }
    ntokens = Tokenizer_split(tkz, str, sep);
    if(verbose) {
        printf("ntokens=%d\n", ntokens);
        fflush(stdout);
    }
    assert(ntokens == Tokenizer_ntokens(tkz));
    assert(ntokens == 5);
    assert(strcmp(Tokenizer_token(tkz, 0), " now") == 0);
    assert(strcmp(Tokenizer_token(tkz, 1), "   \t") == 0);
    assert(strcmp(Tokenizer_token(tkz, 2), " \t is       ") == 0);
    assert(strcmp(Tokenizer_token(tkz, 3), "the  ") == 0);
    assert(strcmp(Tokenizer_token(tkz, 4), "time  \n") == 0);

    if(verbose) {
        for(i = 0; i < ntokens; ++i)
            printf("%4d \"%s\"\n", i, Tokenizer_token(tkz, i));
    }

    ntokens = Tokenizer_strip(tkz, " \t\n");
    assert(ntokens == 4);
    assert(strcmp(Tokenizer_token(tkz, 0), "now") == 0);
    assert(strcmp(Tokenizer_token(tkz, 1), "is") == 0);
    assert(strcmp(Tokenizer_token(tkz, 2), "the") == 0);
    assert(strcmp(Tokenizer_token(tkz, 3), "time") == 0);

    assert(Tokenizer_find(tkz, "now") == 0);
    assert(Tokenizer_find(tkz, "is") == 1);
    assert(Tokenizer_find(tkz, "the") == 2);
    assert(Tokenizer_find(tkz, "time") == 3);
    assert(Tokenizer_find(tkz, "not there") == ntokens);

    if(verbose) {
        printf("after stripping extraneous chars, ntokens is %d\n", ntokens);
        for(i = 0; i < ntokens; ++i)
            printf("%4d \"%s\"\n", i, Tokenizer_token(tkz, i));
        printf("Tokenizer_print:\n");
        Tokenizer_print(tkz, stdout);
    }

    strcpy(str, "afasf");
    ntokens = Tokenizer_split(tkz, str, ":");
    assert(ntokens == 1);

    strcpy(str, "");
    ntokens = Tokenizer_split(tkz, str, ":");
    assert(ntokens == 0);

    Tokenizer_free(tkz);

    unitTstResult("Tokenizer", "OK");

    return 0;
}
