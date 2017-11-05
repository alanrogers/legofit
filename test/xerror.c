/**
 * @file xerror.c
 * @author Alan R. Rogers
 * @brief Test error.c.
 * @copyright Copyright (c) 2017, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "error.h"
#include "misc.h"
#include <stdio.h>
#include <strings.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif
 
int main(void) {

    char buff[100];

    assert(0 == mystrerror_r(NO_ANCESTRAL_ALLELE, buff, sizeof buff));
    assert(0 == strcmp(buff, "No ancestral allele"));

    assert(0 == mystrerror_r(REF_ALT_MISMATCH, buff, sizeof buff));
    assert(0 == strcmp(buff, "Inconsistent REF and ALT alleles"));
    
    assert(0 == mystrerror_r(BUFFER_OVERFLOW, buff, sizeof buff));
    assert(0 == strcmp(buff, "Buffer overflow"));
    
    assert(0 == mystrerror_r(BAD_RAF_INPUT, buff, sizeof buff));
    assert(0 == strcmp(buff, "Bad .raf input file"));
    
    assert(0 == mystrerror_r(BAD_SORT, buff, sizeof buff));
    assert(0 == strcmp(buff, "Incorrect sort"));
    
    unitTstResult("mystrerror_r", "OK");
    
    return 0;
}
