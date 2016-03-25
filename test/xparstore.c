/**
 * @file xparstore.c
 * @author Alan R. Rogers
 * @brief Test parstore.c.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "parstore.h"
#include <stdio.h>
#include <assert.h>

#ifdef NDEBUG
#  error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(void) {
    ParStore *ps = ParStore_new();
    assert(ParStore_nPar(ps) == 0);

    return 0;
}
