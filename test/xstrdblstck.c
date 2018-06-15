/**
 * @file xstrdblstck.c
 * @author Daniel R. Tabin
 * @brief Unit tests for clic
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "strdblstck.h"
#include "misc.h"

int main(int argc, char **argv){
    StrDblStack* a = NULL;
    StrDblStack* b = NULL;
    StrDblStack* c = NULL;
    StrDblStack* d = NULL;
    StrDblStack* e = NULL;
    StrDblStack* f = NULL;

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

    assert(StrDblStack_compare(a,b) == 0);
    assert(StrDblStack_length(a) == 0);

    a = StrDblStack_push(a, "test", 1.0);

    if(verbose){
      fprintf(stderr,"Stack a:\n");
      StrDblStack_print(a,stderr);
      fprintf(stderr,"Stack b:\n");
      StrDblStack_print(b,stderr);
    }

    assert(StrDblStack_compare(a,b) != 0);
    assert(StrDblStack_length(a) == 1);

    a = StrDblStack_pop(a, &temp);
    a = StrDblStack_pop(a, &temp);

    if(verbose){
      fprintf(stderr,"Stack a:\n");
      StrDblStack_print(a,stderr);
      fprintf(stderr,"Stack b:\n");
      StrDblStack_print(b,stderr);
    }

    assert(StrDblStack_compare(a,b) == 0);
    assert(StrDblStack_length(a) == 0);

    c = parseLegofit_CLIC("s1boot0.legofit");
    d = parseLegofit_CLIC("s1boot0.legofit");

    assert(StrDblStack_compare(c,d) == 0);

    if(verbose){
      fprintf(stderr,"Stack c:\n");
      StrDblStack_print(c,stderr);
      fprintf(stderr,"Stack d:\n");
      StrDblStack_print(d,stderr);
    }

    c = StrDblStack_push(c, "test", 1);

    if(verbose){
      fprintf(stderr,"Stack c:\n");
      StrDblStack_print(c,stderr);
      fprintf(stderr,"Stack d:\n");
      StrDblStack_print(d,stderr);
    }

    assert(StrDblStack_compare(c,d) != 0);


    c = StrDblStack_pop(c, &temp);
    d = StrDblStack_pop(d, &temp);

    assert(StrDblStack_compare(c,d) != 0);

    if(verbose){
      fprintf(stderr,"Stack c:\n");
      StrDblStack_print(c,stderr);
      fprintf(stderr,"Stack d:\n");
      StrDblStack_print(d,stderr);
    }

    e = parseSitPat("s1boot0.legofit");
    f = parseSitPat("s1boot0.legofit");

    assert(StrDblStack_compare(e,f) == 0);

    if(verbose){
      fprintf(stderr,"Stack e:\n");
      StrDblStack_print(e,stderr);
      fprintf(stderr,"Stack f:\n");
      StrDblStack_print(f,stderr);
    }

    e = StrDblStack_push(e, "test", 1);

    if(verbose){
      fprintf(stderr,"Stack e:\n");
      StrDblStack_print(e,stderr);
      fprintf(stderr,"Stack f:\n");
      StrDblStack_print(f,stderr);
    }

    assert(StrDblStack_compare(e,f) != 0);


    e = StrDblStack_pop(e, &temp);
    f = StrDblStack_pop(f, &temp);

    assert(StrDblStack_compare(e,f) != 0);

    if(verbose){
      fprintf(stderr,"Stack e:\n");
      StrDblStack_print(e,stderr);
      fprintf(stderr,"Stack f:\n");
      StrDblStack_print(f,stderr);
    }

    printf("All tests completed\n");

    return 0;
}
