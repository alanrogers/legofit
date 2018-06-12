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

    StrDbl* temp = NULL;

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
    if(verbose){
      f = fopen("xstrdblstck_output.txt", "w");
      FILE* f;
    }

    StrDblStack_print(a,f);

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

    a = StrDblStack_pop(a, temp);
    a = StrDblStack_pop(a, temp);

    if(verbose){
      fprintf(stderr,"Stack a:\n");
      StrDblStack_print(a,stderr);
      fprintf(stderr,"Stack b:\n");
      StrDblStack_print(b,stderr);
    }

    assert(StrDblStack_compare(a,b) == 0);
    assert(StrDblStack_length(a) == 0);

    c = parseLegofit("s1boot.legofit");
    d = parseLegofit("s1boot.legofit");

    assert(StrDblStack_compare(c,d) == 0);

    if(verbose){
      fprintf(stderr,"Stack c:\n");
      StrDblStack_print(c,stderr);
      fprintf(stderr,"Stack d:\n");
      StrDblStack_print(d,stderr);
    }

    c = StrDblStack_push(c, buff, 1);

    if(verbose){
      fprintf(stderr,"Stack c:\n");
      StrDblStack_print(c,stderr);
      fprintf(stderr,"Stack d:\n");
      StrDblStack_print(d,stderr);
    }

    assert(StrDblStack_compare(c,d) != 0);


    c = StrDblStack_pop(c, temp);
    d = StrDblStack_pop(d, temp);

    assert(StrDblStack_compare(c,d) != 0);

    if(verbose){
      fprintf(stderr,"Stack c:\n");
      StrDblStack_print(c,stderr);
      fprintf(stderr,"Stack d:\n");
      StrDblStack_print(d,stderr);

      fclose(f);
    }

    printf("All tests completed\n", );

    return 0;
}
