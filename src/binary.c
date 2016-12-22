/**
@file binary.c
@brief Functions for fiddling with bits.

@copyright Copyright (c) 2016, Alan R. Rogers 
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "binary.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

/// Return x after reversing the order of the bits.
tipId_t reverseBits(tipId_t x) {
    if(sizeof(x)==4)
        return rev32(x);
    assert(sizeof(x) == 8);
    return rev64(x);
}

/// Reverse bits in a 32-bit integer.  p 129 of Hacker's Delight, 2nd
/// edn.
uint32_t rev32(uint32_t x) {
    x = (x & 0x55555555) <<  1 | (x & 0xAAAAAAAA) >>  1;
    x = (x & 0x33333333) <<  2 | (x & 0xCCCCCCCC) >>  2;
    x = (x & 0x0F0F0F0F) <<  4 | (x & 0xF0F0F0F0) >>  4;
    x = (x & 0x00FF00FF) <<  8 | (x & 0xFF00FF00) >>  8;
    x = (x & 0x0000FFFF) << 16 | (x & 0xFFFF0000) >> 16;
    return x;
}

/// Reverse bits in a 64-bit integer.  This just extends the pattern
/// in rev32.
uint64_t rev64(uint64_t x) {
    x = (x & 0x5555555555555555) <<  1
        | (x & 0xAAAAAAAAAAAAAAAA) >>  1;
    x = (x & 0x3333333333333333) <<  2
        | (x & 0xCCCCCCCCCCCCCCCC) >>  2;
    x = (x & 0x0F0F0F0F0F0F0F0F) <<  4
        | (x & 0xF0F0F0F0F0F0F0F0) >>  4;
    x = (x & 0x00FF00FF00FF00FF) <<  8
        | (x & 0xFF00FF00FF00FF00) >>  8;
    x = (x & 0x0000FFFF0000FFFF) << 16
        | (x & 0xFFFF0000FFFF0000) >> 16;
    x = (x & 0x00000000FFFFFFFF) << 32
        | (x & 0xFFFFFFFF00000000) >> 32;
    return x;
}

/// Print the bits in an object of size "size", pointed to by
/// "ptr". Assumes little endian   
void printBits(size_t size, void const * const ptr, FILE *fp) {
    unsigned char const * const b = (unsigned char const * const) ptr;
    unsigned char byte;
    size_t i, j;

    i = size;
    while(i-- > 0) {
        j=8;
        while(j-- > 0) {
            byte = b[i] & (1<<j);
            byte >>= j;
            fprintf(fp, "%u", (unsigned) byte);
        }
    }
    putc('\n', fp);
}

/**
 * Examine the bits in x. If bit i equals 1, then add i to array bit.
 * Return the number of nonzero bits. maxbits is the allocated length
 * of array bit. It is an error if the number of nonzero bits exceeds
 * maxbits.
 */
int getBits(tipId_t x, int maxbits, int bit[maxbits]) {
    size_t len = 8*sizeof(x);
    size_t i;
    int nbits=0;

    for(i=0; i < len; ++i) {
        if( (x>>i) & 1) {
            if(nbits == maxbits) {
                fprintf(stderr,"%s:%d: Too many nonzero bits in %lu."
                        " Max is %d\n",
                        __FILE__, __LINE__, (unsigned long) x, maxbits);
                exit(1);
            }
            bit[nbits++] = i;
        }
    }
    return nbits;
}

/// Print the (0-based) indices of the bits that are "on".
/// @param [in] ptr points to the object whose bits will be examined
/// @param [in] size the size of the object in bytes
void printWhichBits(size_t const size, void const * const ptr) {
    unsigned char const * const b = (unsigned char const * const) ptr;
    unsigned char byte;
    size_t i, j;
    int needComma=0;

    for(i=0; i < size; ++i) {
        for(j=0; j < 8; ++j) {
            byte = b[i] & (1<<j);
            if(byte) {
                if(needComma)
                    putchar(',');
                else
                    needComma = 1;
                printf("%zu", 8*i + j);
            }
        }
    }
}

/// Count the number of 1 bits in x. From p 85 of Hacker's Delight,
/// 2nd Edn, by Henry S. Warren, Jr. This algorithm is fast when the
/// number of 1 bits is small, and it makes no assumption about the
/// number of bits in x.
/// @return the number of bits.
int num1bits(tipId_t x) {
    int n=0;

    while(x != 0) {
        ++n;
        x = x & (x-1);
    }
    return n;
}
