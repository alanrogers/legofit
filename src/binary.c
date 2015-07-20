#include "binary.h"
#include <stdio.h>
#include <stdlib.h>

/* assumes little endian */
void printBits(size_t size, void const * const ptr) {
    unsigned char const * const b = (unsigned char const * const) ptr;
    unsigned char byte;
    size_t i, j;

    i = size;
    while(i-- > 0) {
        j=8;
        while(j-- > 0) {
            byte = b[i] & (1<<j);
            byte >>= j;
            printf("%u", (unsigned) byte);
        }
    }
    putchar('\n');
}

/*
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

/// Count the number of 1 bits. From p 85 of Hacker's Delight, 2nd
/// Edn, by Henry S. Warren, Jr. This algorithm is fast when the
/// number of 1 bits is small, and it makes no assumption about the 
/// number of bits in x.
int num1bits(tipId_t x) {
    int n=0;

    while(x != 0) {
        ++n;
        x = x & (x-1);
    }
    return n;
}
