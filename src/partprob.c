#include "partprob.h"
#include <math.h>
#include <assert.h>
/**
 * @file partprob.c
 * @author Alan R. Rogers
 * @brief Probability of a set partition
 *
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

/** 
 * Log of constant in coalescent probability from theorem 1.5, p. 11, Durrett,
 * Richard. 2008. Probability Models for DNA Sequence Evolution.
 */
long double lnCoalConst(unsigned n, unsigned k) {
    assert(n >= k);
    assert(k > 0);
    return lgammal(k+1)
        - lgammal(n+1)
        + lgammal(n-k+1)
        + lgammal(k)
        - lgammal(n);
}

/**
 * Probability of a set partition. There are n descendants in some
 * recent epoch and k < n in some earlier epoch. In that earlier
 * epoch, the i'th ancestor had y[i] descendants. The function returns
 * the probability of a partition that satisfies this condition, given
 * n and k. There may be several such partitions. lnconst should be
 * calculated using function lnCoalConst. See theorem 1.5, p. 11,
 * Durrett, Richard. 2008. Probability Models for DNA Sequence
 * Evolution.
 */
double probPartition(unsigned k, unsigned y[k], long double lnconst) {
    assert(k > 0);
    long double x = 0.0L;
    for(unsigned i=0; i<k; ++i)
        x += lgammal(y[i] + 1);
    return (double) expl(lnconst + x);
}

