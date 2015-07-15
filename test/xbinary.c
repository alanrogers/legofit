/**
 * @file xbinary.c
 * @brief Test file binary.c
 *
 * @copyright Copyright (c) 2014, Alan R. Rogers 
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "binary.h"
#include <stdio.h>
#include <limits.h>

int main(int argv, char* argc[])
{
    int i = 2;
    unsigned ui = UINT_MAX;
    float f = 23.45f;
    printBits(sizeof(i), &i);
    printWhichBits(sizeof(i), &i);
    putchar('\n');
    printBits(sizeof(ui), &ui);
    printWhichBits(sizeof(ui), &ui);
    putchar('\n');
    printBits(sizeof(f), &f);
    printWhichBits(sizeof(f), &f);
    putchar('\n');
    return 0;
}
