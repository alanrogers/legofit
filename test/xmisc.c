/**
 * @file xmisc.c
 * @author Alan R. Rogers
 * @brief Test parstore.c.
 * @copyright Copyright (c) 2017, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#ifdef NDEBUG
#  error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {
	int verbose=0;

	switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0) {
            fprintf(stderr, "usage: xmisc [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xmisc [-v]\n");
        exit(EXIT_FAILURE);
    }

    int i;
    char buff[200];
    char *tok, *next = buff;
    strcpy(buff, "  tok0	tok2  	  tok2  ");

    for(i=0; i<3; ++i) {
        tok = nextWhitesepToken(&next);
        assert(tok);
        if(verbose)
            printf("token %d: \"%s\"\n", i, tok);
    }

    tok = nextWhitesepToken(&next);
    assert(tok==NULL);

    unitTstResult("nextWhitesepToken", "OK");

    FILE *fp = efopen("xmisc.c", "r");
    fclose(fp);

# if 0
    // This should abort.
    fp = efopen("NotThere", "r");
#endif

    unitTstResult("efopen", "OK");

    // Test tokenize
    int n, dim=10;
    char s[100];
    char *tokens[dim];

    strcpy(s, "aaa,bbb,c");
    n = tokenize(dim, tokens, s, ",");
    assert(n==3);
    assert(strcmp(tokens[0], "aaa") == 0);
    assert(strcmp(tokens[1], "bbb") == 0);
    assert(strcmp(tokens[2], "c") == 0);
    if(verbose) {
        for(i=0; i<n; ++i)
            printf("token[%d] = %s\n", i, tokens[i]);
    }

    strcpy(s, "aaa");
    n = tokenize(dim, tokens, s, ",");
    assert(n==1);
    assert(strcmp(tokens[0], "aaa") == 0);
    if(verbose) {
        for(i=0; i<n; ++i)
            printf("token[%d] = %s\n", i, tokens[i]);
    }

    unitTstResult("tokenize", "OK");

    // Test strReplaceChr
    strcpy(s, "a-b-c");
    strReplaceChr(s, '-', '.');
    assert(0==strcmp(s, "a.b.c"));

    unitTstResult("strReplaceChr", "OK");

    // test parseDbl
    double x;
    strcpy(buff, " 123.4 ");
    errno=0;
    x = parseDbl(buff);
    assert(errno==0);
    assert(Dbl_near(x, 123.4));

    strcpy(buff, " 123.4e1 ");
    errno=0;
    x = parseDbl(buff);
    assert(errno==0);
    assert(Dbl_near(x, 123.4e1));

    strcpy(buff, " 123.4e9999 ");
    errno=0;
    x = parseDbl(buff);
    assert(errno==ERANGE);
    assert(x==0.0);

    strcpy(buff, " 123.4xxxx ");
    errno=0;
    x = parseDbl(buff);
    assert(errno==EINVAL);
    assert(x==0.0);

    strcpy(buff, " xxxx ");
    errno=0;
    x = parseDbl(buff);
    assert(errno==EINVAL);
    assert(x==0.0);
    
    unitTstResult("parseDbl", "OK");

    return 0;
}
