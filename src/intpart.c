/**
   @file intpart.c
   @brief Partitions of an integer into a given number of summands.

   @copyright Copyright (c) 2020, Alan R. Rogers
   <rogers@anthro.utah.edu>. This file is released under the Internet
   Systems Consortium License, which can be found in file "LICENSE".
*/
#include "intpart.h"
#include "misc.h"
#include "u64u64map.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

// So we don't have to calculate the same value more than once.
static pthread_mutex_t map_lock = PTHREAD_MUTEX_INITIALIZER;
static U64U64Map *map=NULL;

static uint64_t numIntPart_r(int32_t n, int32_t k);

/// Number of ways to partition a positive integer n into k parts. The
/// public-facing function just locks the map and calls the static
/// recursive function, numIntPart_r.
uint64_t numIntPart(int32_t n, int32_t k) {
    int lockstat = pthread_mutex_lock(&map_lock);
    if(lockstat)
        ERR(lockstat, "lock");

    if(map == NULL) {
        // allocate hash table on first call
        map = U64U64Map_new(512);
        CHECKMEM(map);
    }

    uint64_t value = numIntPart_r(n, k);

    lockstat = pthread_mutex_unlock(&map_lock);
    if(lockstat)
        ERR(lockstat, "unlock");

    return value;
}

/// Called by numIntPart. This assumes that the map is locked and does
/// the real work of calculating the number of ways to partition a
/// positive integer n into k parts.
static uint64_t numIntPart_r(int32_t n, int32_t k) {
    uint64_t key, value;
    int status;

    if(n==0 && k==0)
        return 1ULL;
    
    if(n<=0 || k<=0)
        return 0ULL;

    // form 64-bit key from two 32-bit arguments
    key = (uint32_t) n;
    key <<= 32;
    key |= (uint32_t) k;

    value = U64U64Map_get(map, key, &status);
    if(status == 0)
        return value;

    value = numIntPart_r(n-k, k) + numIntPart_r(n-1, k-1);
    status = U64U64Map_insert(map, key, value);
    if(status) {
        fprintf(stderr,"%s:%d: inserted duplicated value\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    return value;
}

/// Partition a positive integer n into a sum of k positive integers.
/// Algorithm H, p 392 of Knuth, Donald E. 2011. The art of computer
/// programming, volume 4A.
int traverseIntPartitions(int n, int k,
                          int (*visit)(int kk, int a[kk], void *data),
                          void *data) {
    if(k < 2) {
        fprintf(stderr,"%s:%d: k (%d) must be > 2\n",
                __FILE__,__LINE__, k);
        exit(EXIT_FAILURE);
    }
    if(k > n) {
        fprintf(stderr,"%s:%d: k (%d) must be <= n (%d)\n",
                __FILE__,__LINE__, k, n);
        exit(EXIT_FAILURE);
    }
    int j;
    int a[k+1];
    a[0] = n - k + 1;
    for(j=1; j < k; ++j)
        a[j] = 1;
    a[k] = -1;

    while(1) {
        int status = visit(k, a, data);
        if(status)
            return status;
        if(a[1] < a[0] - 1) {
            a[0] -= 1;
            a[1] += 1;
            continue;
        }

        j = 2;
        int s = a[0] + a[1] - 1;
        while(a[j] >= a[0] - 1) {
            s += a[j];
            j += 1;
        }
        if(j+1>k)
            break;
        int x = a[j] + 1;
        a[j] = x;
        j -= 1;
        while(j>0) {
            a[j] = x;
            s -= x;
            j -= 1;
        }
        a[0] = s;
    }
    return 0;
}

/// Free static variable.
void numIntPart_free(void) {
    int status = pthread_mutex_lock(&map_lock);
    if(status)
        ERR(status, "lock");

    if(map) {
        U64U64Map_free(map);
        map = NULL;
    }

    status = pthread_mutex_unlock(&map_lock);
    if(status)
        ERR(status, "unlock");
}
