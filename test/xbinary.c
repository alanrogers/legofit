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
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

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

    ui = UINT_MAX;
    printf("%12s: ", "UINT_MAX");
    printBits(sizeof(ui), &ui);
    ui >>= 2;
    printf("%12s: ", "UINT_MAX>>2");
    printBits(sizeof(ui), &ui);

    i = -1;
    printf("%12s: ", "-1");
    printBits(sizeof(i), &i);
    i >>= 2;
    printf("%12s: ", "-1 >>2");
    printBits(sizeof(i), &i);

    ui = -1;
    printf("%12s: ", "-1u");
    printBits(sizeof(ui), &ui);
    ui >>= 2;
    printf("%12s: ", "-1u >>2");
    printBits(sizeof(ui), &ui);

    assert(!isPow2(0UL));
    assert(isPow2(1UL));
    assert(isPow2(2UL));
    assert(!isPow2(3UL));
    assert(isPow2(4UL));
    assert(!isPow2(5UL));
    assert(!isPow2(6UL));
    assert(!isPow2(7UL));
    assert(isPow2(8UL));

    assert(num1bits(0UL)==0);
    assert(num1bits(1UL)==1);
    assert(num1bits(2UL)==1);
    assert(num1bits(3UL)==2);
    assert(num1bits(4UL)==1);
    assert(num1bits(5UL)==2);
    assert(num1bits(6UL)==2);
    assert(num1bits(7UL)==3);
    assert(num1bits(8UL)==1);
    assert(num1bits((1UL << 1) | (1UL<<3) | (1UL<<5) | (1UL<<9))==4);

    return 0;
}
