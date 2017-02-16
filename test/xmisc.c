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

    return 0;
}
