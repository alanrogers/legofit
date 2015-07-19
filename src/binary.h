#ifndef ARR_BINARY_H
#define ARR_BINARY_H

#include <stddef.h>

void printBits(size_t size, void const * const ptr);
void printWhichBits(size_t size, void const * const ptr);
int getBits(unsigned long x, int maxbits, int *bit);
int num1bits(unsigned long x);
static inline int isPow2(unsigned long x);

static inline int isPow2(unsigned long x) {
    return (x>0 && !(x & (x-1)));
}

#endif
