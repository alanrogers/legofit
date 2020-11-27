#ifndef ARR_BINARY_H
#  define ARR_BINARY_H

#  include "typedefs.h"
#  include <stddef.h>
#  include <stdint.h>
#  include <stdio.h>

void        printBits(size_t size, void const *const ptr, FILE *fp);
void        printWhichBits(size_t size, void const *const ptr);
int         getBits(tipId_t x, int maxbits, int *bit);
int         num1bits(tipId_t x);
static inline tipId_t reverseBits(tipId_t x);
uint32_t    rev32(uint32_t x);
uint64_t    rev64(uint64_t x);
#ifdef __GNUC__
uint32_t    uint32Hash(uint32_t key );
uint32_t    uint64Hash(uint64_t key);
#else
uint32_t    uint32Hash(uint32_t key ) __attribute__((no_sanitize("integer")));
uint32_t    uint64Hash(uint64_t key) __attribute__((no_sanitize("integer")));
#endif
int         nlz(tipId_t x);
int         nlz64(uint64_t x);
int         no_shared_bits(int n, tipId_t *tid);
static inline int isPow2(tipId_t x);
tipId_t     next_power_of_2(tipId_t x);
tipId_t     low_bits_on(unsigned n);

static inline int isPow2(tipId_t x) {
    return (x > 0 && !(x & (x - 1)));
}

/// Return x after reversing the order of the bits.
static inline tipId_t reverseBits(tipId_t x) {
#if TIPID_SIZE==32
    return rev32(x);
#elif TIPID_SIZE==64    
    return rev64(x);
#else
#   error "unsupported tipId_t size"
#endif    
}

#endif
