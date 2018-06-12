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
    char buff[200] = {'t', 'e', 's', 't', '\0'};

    f = fopen("output.txt", "w");

    StrDblStack_print(a,f);

    assert(StrDblStack_compare(a,b) == 0);
    assert(StrDblStack_length(a) == 0);

    StrDblStack_push(a, buff, 1);

    assert(StrDblStack_compare(a,b) != 0);
    assert(StrDblStack_length(a) == 1);

    StrDblStack_pop(a, temp);
    StrDblStack_pop(a, temp);

    assert(StrDblStack_compare(a,b) == 0);
    assert(StrDblStack_length(a) == 0);

    c = parseLegofit("s1boot.legofit");
    d = parseLegofit("s1boot.legofit");

    assert(StrDblStack_compare(c,d) == 0);

    StrDblStack_push(c, buff, 1);

    assert(StrDblStack_compare(c,d) != 0);

    StrDblStack_print(c,f);

    StrDblStack_pop(c, temp);
    StrDblStack_pop(d, temp);

    assert(StrDblStack_compare(c,d) != 0);

    return 0;
}
