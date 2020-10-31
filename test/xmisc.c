/**
 * @file xmisc.c
 * @author Alan R. Rogers
 * @brief Test parstore.c.
 * @copyright Copyright (c) 2017, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "misc.h"
#include "error.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <unistd.h>

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
    fclose(fp);

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

    strcpy(buff, "abcdefghijklmn");
    assert(0 == strcmp("lmn", strltrunc(buff, 3)));

    strcpy(buff, "abc");
    assert(0 == strcmp("abc", strltrunc(buff, 3)));

    strcpy(buff, "abc");
    assert(0 == strcmp("abc", strltrunc(buff, 4)));

    unitTstResult("strltrunc", "OK");

    unsigned y[] = {0, 0, 1, 0, 0, 0, 2, 0, 0, 3, 0, 0};
    dim = (int) (sizeof(y)/sizeof(y[0]));
    dim = removeZeroes(dim, y);
    assert(dim == 3);
    for(i=0; i < dim; ++i) {
        if(verbose)
            printf("%u ", y[i]);
        assert(y[i] == i+1u);
    }
    if(verbose)
        putchar('\n');
    unitTstResult("removeZeroes", "OK");

    fp = fopen("xmisc.tmp", "w");
    assert(fp);
    fputs("123456789", fp);
    fclose(fp);

    const char *tmpfile = "xmisc.tmp";
    int status;
    fp=fopen(tmpfile, "r");
    status = readline(5, buff, fp);
    assert(status == BUFFER_OVERFLOW);
    rewind(fp);

    status = readline(9, buff, fp);
    assert(status == BUFFER_OVERFLOW);
    rewind(fp);

    status = readline(10, buff, fp);
    assert(status == 0);

    status = readline(10, buff, fp);
    assert(status == EOF);
    
    unitTstResult("readline", "OK");

    unlink(tmpfile);

    const char *str = mybasename("/a/b/cde");
    assert(0 == strcmp(str, "cde"));
    str = mybasename("a/b/cde/");
    assert(strlen(str) == 0);

    unitTstResult("mybasename", "OK");

    char dst[5];
    const char *src = "abcd";
    status = strnncopy(sizeof(dst), dst, 3, src);
    assert(status==0);
    assert(0 == strcmp(dst, "abc"));
    assert(strlen(dst) == 3);

    src = "abcde";
    status = strnncopy(sizeof(dst), dst, 5, src);
    assert(status==BUFFER_OVERFLOW);
    assert(strlen(dst) < sizeof(dst));

    src = "abcdef";
    status = strnncopy(sizeof(dst), dst, 6, src);
    assert(status==BUFFER_OVERFLOW);
    assert(strlen(dst) < sizeof(dst));

    src = "abcdef";
    status = strnncopy(sizeof(dst), dst, 0, src);
    assert(status==0);
    assert(strlen(dst) == 0);

    unitTstResult("strnncopy", "OK");

    // test strchrcnt
    assert(1 == strchrcnt("asd", '\0'));
    assert(1 == strchrcnt("asd", 'a'));
    assert(0 == strchrcnt("asd", 'e'));
    assert(3 == strchrcnt("a,s,,d", ','));
    assert(4 == strchrcnt(",,,,", ','));
    
    unitTstResult("strchrcnt", "OK");

    // test collapse_whitespace
    strcpy(buff, " aa   bcde                f ");
    collapse_whitespace(buff);
    assert(0 == strcmp(buff, " aa bcde f "));

    strcpy(buff, "aa");
    collapse_whitespace(buff);
    assert(0 == strcmp(buff, "aa"));

    strcpy(buff, "");
    collapse_whitespace(buff);
    assert(0 == strcmp(buff, ""));

    strcpy(buff, "          ");
    collapse_whitespace(buff);
    assert(0 == strcmp(buff, " "));

    unitTstResult("collapse_whitespace", "OK");
        
    return 0;
}
