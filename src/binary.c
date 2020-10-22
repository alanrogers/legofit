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
#include <assert.h>

static int nlz32(uint32_t x);

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

/// Number of leading zeroes in 32-bit unsigned integer.  From p 99 of
/// Hacker's Delight, 2nd edition.
static int nlz32(uint32_t x) {
    if(x == 0)
        return 32;
    int n = 1;
    if((x >> 16) == 0) {n += 16; x <<= 16;}
    if((x >> 24) == 0) {n += 8; x <<= 8;}
    if((x >> 28) == 0) {n += 4; x <<= 4;}
    if((x >> 30) == 0) {n += 2; x <<= 2;}
    n -= (x >> 31);
    return n;
}

/// Number of leading zeroes in 64-bit unsigned integer.  64-bit
/// version of algorithm on p 99 of Hacker's Delight, 2nd edition.
int nlz64(uint64_t x) {
    if(x == 0)
        return 64;
    int n = 1;
    if((x>>32) == 0) {n += 32; x <<= 32;};
    if((x>>48) == 0) {n += 16; x <<= 16;};
    if((x>>56) == 0) {n += 8; x <<= 8;};
    if((x>>60) == 0) {n += 4; x <<= 4;};
    if((x>>62) == 0) {n += 2; x <<= 2;};
    n -= (x >> 63);
    return n;
}

/// Number of leading zeroes in tipId_t variable.
int nlz(tipId_t x) {
#if TIPID_SIZE==32    
    return nlz32(x);
#elif TIPID_SIZE==64    
    return nlz64(x);
#else
#   error "unsupported tipId_t size"
#endif    
}

/// Print the bits in an object of size "size", pointed to by
/// "ptr". Assumes little endian   
void printBits(size_t size, void const * const ptr, FILE *fp) {
    unsigned char const * const b = (unsigned char const * const) ptr;
    unsigned char byte;
    size_t i, j;

    i = size;
    while(i > 0) {
        i -= 1;
        j=8;
        while(j > 0) {
            j -= 1;
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

/// Round up to the next largest power of 2.
tipId_t next_power_of_2(tipId_t x) {
    if(x == 0)
        return 1;
    if( (x & (x-1)) == 0 )
        return x; // x is already a power of 2

    // Turn off all bits except the highest.
    do {
        x &= x-1;
    }while(x & (x-1));

    // Shift the highest bit one position to the left.
    return x << 1;
}

/// Hash function for a 32-bit integer. From Thomas Wang's 1997
/// article: 
/// https://gist.github.com/badboy/6267743
/// This generates integer overflows, presumably by design.
#ifdef __clang__
uint32_t uint32Hash( uint32_t key ) __attribute__((no_sanitize("integer")))
#else
    uint32_t uint32Hash( uint32_t key )
#endif
{    
    key = (key+0x7ed55d16) + (key<<12);
    key = (key^0xc761c23c) ^ (key>>19);
    key = (key+0x165667b1) + (key<<5);
    key = (key+0xd3a2646c) ^ (key<<9);
    key = (key+0xfd7046c5) + (key<<3);
    key = (key^0xb55a4f09) ^ (key>>16);
    return key;
}

/// Hash function for a 64-bit integer.
#ifdef __clang__
uint32_t uint64Hash(uint64_t key) __attribute__((no_sanitize("integer")))
#else
uint32_t uint64Hash(uint64_t key)
#endif
{    
    key = (~key) + (key << 18);
    key = key ^ (key >> 31);
    key = key * 21;
    key = key ^ (key >> 11);
    key = key + (key << 6);
    key = key ^ (key >> 22);
    return (uint32_t) key;
}

/// Return 1 if the values in array share no bits; 0 otherwise.
/// Returns 1 if n==0.
int no_shared_bits(int n, tipId_t *tid) {
    tipId_t u = 0; // union of prior tipId_t values
    for(int i=0; i < n; ++i) {
        if(u & tid[i])
            return 0;
        u |= tid[i];
    }
    return 1;
}

