#ifndef ARR_BINARY_H
#define ARR_BINARY_H

#include "typedefs.h"
#include <stddef.h>

void printBits(size_t size, void const * const ptr);
void printWhichBits(size_t size, void const * const ptr);
int getBits(tipId_t x, int maxbits, int *bit);
int num1bits(tipId_t x);
static inline int isPow2(tipId_t x);

static inline int isPow2(tipId_t x) {
    return (x>0 && !(x & (x-1)));
}

#endif
