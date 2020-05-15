#ifndef ARR_BINARY_H
#  define ARR_BINARY_H

#  include "typedefs.h"
#  include <stddef.h>
#  include <stdio.h>

void        printBits(size_t size, void const *const ptr, FILE *fp);
void        printWhichBits(size_t size, void const *const ptr);
int         getBits(tipId_t x, int maxbits, int *bit);
int         num1bits(tipId_t x);
tipId_t     reverseBits(tipId_t x);
uint32_t    rev32(uint32_t x);
uint64_t    rev64(uint64_t x);
uint32_t    uint32Hash(uint32_t key ) __attribute__((no_sanitize("integer")));
uint32_t    uint64Hash(uint64_t key) __attribute__((no_sanitize("integer")));
int         nlz(tipId_t x);
int         nlz64(uint64_t x);
static inline int isPow2(tipId_t x);

static inline int isPow2(tipId_t x) {
    return (x > 0 && !(x & (x - 1)));
}

#endif
