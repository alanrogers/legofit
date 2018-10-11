/**
 * @file xbinary.c
 * @brief Test file binary.c
 *
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "binary.h"
#include "misc.h"
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {

    int         i, verbose = 0;

    switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0) {
            fprintf(stderr, "usage: xbinary [-v]\n");
            exit(1);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xbinary [-v]\n");
        exit(1);
    }

    i = 2;
    unsigned ui = UINT_MAX;
    float f = 23.45f;
    if(verbose) {
        printBits(sizeof(i), &i, stdout);
        printWhichBits(sizeof(i), &i);
        putchar('\n');
        printBits(sizeof(ui), &ui, stdout);
        printWhichBits(sizeof(ui), &ui);
        putchar('\n');
        printBits(sizeof(f), &f, stdout);
        printWhichBits(sizeof(f), &f);
        putchar('\n');
    }

    ui = UINT_MAX;
    if(verbose) {
        printf("%12s: ", "UINT_MAX");
        printBits(sizeof(ui), &ui, stdout);
    }
    ui >>= 2;
    if(verbose) {
        printf("%12s: ", "UINT_MAX>>2");
        printBits(sizeof(ui), &ui, stdout);
    }

    i = -1;
    if(verbose) {
        printf("%12s: ", "-1");
        printBits(sizeof(i), &i, stdout);
    }
    i >>= 2;
    if(verbose) {
        printf("%12s: ", "-1 >>2");
        printBits(sizeof(i), &i, stdout);
    }

    ui = -1;
    if(verbose) {
        printf("%12s: ", "-1u");
        printBits(sizeof(ui), &ui, stdout);
    }
    ui >>= 2;
    if(verbose) {
        printf("%12s: ", "-1u >>2");
        printBits(sizeof(ui), &ui, stdout);
    }

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

    uint32_t u32, r32;
    u32 = 1 | (1 << 3) | (1 << 17);
    r32 = rev32(u32);
    assert(u32 == rev32(r32));
    if(verbose) {
        printf("u32:");
        printBits(sizeof(u32), &u32, stdout);
        putchar('\n');
        printf("r32:");
        printBits(sizeof(r32), &r32, stdout);
        putchar('\n');
    }

    uint64_t u64, r64;
    u64 = 1ULL | (1ULL << 3) | (1ULL << 17) | (1ULL << 47);
    r64 = rev64(u64);
    assert(u64 == rev64(r64));
    if(verbose) {
        printf("u64:");
        printBits(sizeof(u64), &u64, stdout);
        putchar('\n');
        printf("r64:");
        printBits(sizeof(r64), &r64, stdout);
        putchar('\n');
    }

    uint32_t key32 = 1234u;
    if(verbose)
        printf("32-bit key %u -> hash %u\n", key32, uint32Hash(key32));

    uint64_t key64 = 1234u;
    if(verbose)
        printf("64-bit key %llu -> hash %u\n", (long long unsigned) key64,
               uint64Hash(key64));

    unitTstResult("binary", "OK");
    return 0;
}
