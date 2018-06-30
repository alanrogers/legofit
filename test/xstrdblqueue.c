/**
 * @file xstrdblqueue.c
 * @author Daniel R. Tabin
 * @brief Unit tests for clic
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "strdblqueue.h"
#include "misc.h"

int main(int argc, char **argv){
    StrDblQueue* a = NULL;
    StrDblQueue* b = NULL;
    StrDblQueue* c = NULL;
    StrDblQueue* d = NULL;
    StrDblQueue* e = NULL;
    StrDblQueue* f = NULL;

    StrDbl temp;

    bool verbose = 0;

	  switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0) {
            fprintf(stderr, "usage: xstrdblstck [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xstrdblstck [-v]\n");
        exit(EXIT_FAILURE);
    }

    assert(StrDblQueue_compare(a,b) == 0);
    assert(StrDblQueue_length(a) == 0);

    a = StrDblQueue_push(a, "test", 1.0);

    if(verbose){
      fprintf(stderr,"Queue a:\n");
      StrDblQueue_print(a,stderr);
      fprintf(stderr,"Queue b:\n");
      StrDblQueue_print(b,stderr);
    }

    assert(StrDblQueue_compare(a,b) != 0);
    assert(StrDblQueue_length(a) == 1);

    a = StrDblQueue_pop(a, &temp);
    a = StrDblQueue_pop(a, &temp);

    if(verbose){
      fprintf(stderr,"Queue a:\n");
      StrDblQueue_print(a,stderr);
      fprintf(stderr,"Queue b:\n");
      StrDblQueue_print(b,stderr);
    }

    assert(StrDblQueue_compare(a,b) == 0);
    assert(StrDblQueue_length(a) == 0);

    c = parseLegofit_CLIC("s1boot0.legofit");
    d = parseLegofit_CLIC("s1boot0.legofit");

    assert(StrDblQueue_compare(c,d) == 0);

    if(verbose){
      fprintf(stderr,"Queue c:\n");
      StrDblQueue_print(c,stderr);
      fprintf(stderr,"Queue d:\n");
      StrDblQueue_print(d,stderr);
    }

    c = StrDblQueue_push(c, "test", 1);

    if(verbose){
      fprintf(stderr,"Queue c:\n");
      StrDblQueue_print(c,stderr);
      fprintf(stderr,"Queue d:\n");
      StrDblQueue_print(d,stderr);
    }

    assert(StrDblQueue_compare(c,d) != 0);


    c = StrDblQueue_pop(c, &temp);
    d = StrDblQueue_pop(d, &temp);

    assert(StrDblQueue_compare(c,d) != 0);

    if(verbose){
      fprintf(stderr,"Queue c:\n");
      StrDblQueue_print(c,stderr);
      fprintf(stderr,"Queue d:\n");
      StrDblQueue_print(d,stderr);
    }

    e = parseSitPat("s1boot0.legofit");
    f = parseSitPat("s1boot0.legofit");

    assert(StrDblQueue_compare(e,f) == 0);

    if(verbose){
      fprintf(stderr,"Queue e:\n");
      StrDblQueue_print(e,stderr);
      fprintf(stderr,"Queue f:\n");
      StrDblQueue_print(f,stderr);
    }

    e = StrDblQueue_push(e, "test", 1);

    if(verbose){
      fprintf(stderr,"Queue e:\n");
      StrDblQueue_print(e,stderr);
      fprintf(stderr,"Queue f:\n");
      StrDblQueue_print(f,stderr);
    }

    assert(StrDblQueue_compare(e,f) != 0);


    e = StrDblQueue_pop(e, &temp);
    f = StrDblQueue_pop(f, &temp);

    assert(StrDblQueue_compare(e,f) != 0);

    if(verbose){
      fprintf(stderr,"Queue e:\n");
      StrDblQueue_print(e,stderr);
      fprintf(stderr,"Queue f:\n");
      StrDblQueue_print(f,stderr);
    }

    printf("All tests completed\n");

    return 0;
}
