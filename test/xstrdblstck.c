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

    FILE* f;

    StrDbl* temp = NULL;

    f = fopen("xstrdblstck_output.txt", "w");

    StrDblStack_print(a,f);

    assert(StrDblStack_compare(a,b) == 0);
    assert(StrDblStack_length(a) == 0);

    a = StrDblStack_push(a, "test.a", 1.0);

    fprintf(stderr,"Stack a:\n");
    StrDblStack_print(a,stderr);
    fprintf(stderr,"Stack b:\n");
    StrDblStack_print(b,stderr);

    assert(StrDblStack_compare(a,b) != 0);
    assert(StrDblStack_length(a) == 1);

    StrDblStack_pop(a, temp);
    StrDblStack_pop(a, temp);

    assert(StrDblStack_compare(a,b) == 0);
    assert(StrDblStack_length(a) == 0);

    c = parseLegofit("s1boot.legofit");
    d = parseLegofit("s1boot.legofit");

    assert(StrDblStack_compare(c,d) == 0);

    c = StrDblStack_push(c, buff, 1);

    assert(StrDblStack_compare(c,d) != 0);

    StrDblStack_print(c,f);

    StrDblStack_pop(c, temp);
    StrDblStack_pop(d, temp);

    assert(StrDblStack_compare(c,d) != 0);

    fclose(f);

    return 0;
}
