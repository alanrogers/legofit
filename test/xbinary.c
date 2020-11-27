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

    // Test nlz
    tipId_t tid;
    if(verbose) {
        for(tid = 0; tid <= 16; ++tid)
            printf("nlz(%#08x) = %d\n", tid, nlz(tid));
    }

    assert( 32 == nlz(0) );
    assert( 31 == nlz(1) );
    assert( 30 == nlz(2) );
    assert( 30 == nlz(3) );
    assert( 29 == nlz(4) );
    assert( 29 == nlz(5) );
    assert( 29 == nlz(6) );
    assert( 28 == nlz(8) );
    assert( 28 == nlz(9) );
    assert( 27 == nlz(16) );
    assert( 27 == nlz(18) );
    assert( 26 == nlz(32) );
    assert( 25 == nlz(64) );
    assert( 24 == nlz(128) );
    assert( 23 == nlz(256) );
    assert(  0 == nlz(0xffffffff) );

    assert( 64 == nlz64(0) );
    assert( 63 == nlz64(1) );
    assert( 62 == nlz64(2) );
    assert( 62 == nlz64(3) );
    assert( 61 == nlz64(4) );
    assert( 61 == nlz64(5) );
    assert( 61 == nlz64(6) );
    assert( 60 == nlz64(8) );
    assert( 60 == nlz64(9) );
    assert( 59 == nlz64(16) );
    assert( 59 == nlz64(18) );
    assert( 58 == nlz64(32) );
    assert( 57 == nlz64(64) );
    assert( 56 == nlz64(128) );
    assert( 55 == nlz64(256) );
    assert( 32 == nlz64(0xffffffffLU) );
    assert( 0 == nlz64(0xffffffffffffffffLU) );
    
    unitTstResult("nlz", "OK");

    assert(1 == next_power_of_2(0));
    assert(1 == next_power_of_2(1));
    assert(2 == next_power_of_2(2));
    assert(4 == next_power_of_2(3));
    assert(4 == next_power_of_2(4));
    assert(8 == next_power_of_2(5));
    assert(8 == next_power_of_2(6));
    assert(8 == next_power_of_2(7));
    assert(8 == next_power_of_2(8));
    assert(16 == next_power_of_2(9));

    unitTstResult("next_power_of_2", "OK");

    unsigned twopwr = 1;
    for(unsigned u=1; u < 10; ++u) {
        twopwr *= 2;
        tipId_t maxtid = low_bits_on(u);
        if(verbose) {
            printf("bits(%u) = ", u);
            printBits(sizeof(tipId_t), &maxtid, stdout);
            putchar('\n');
        }
        assert(maxtid == (1u << u) - 1);
        assert(maxtid == twopwr - 1);
    }
    
    unitTstResult("low_bits_on", "OK");

    unitTstResult("binary", "OK");
    
    return 0;
}
