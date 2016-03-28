/**
 * @file xparkeyval.c
 * @author Alan R. Rogers
 * @brief Test parkeyval.c.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "parkeyval.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char* argv[]){

	int verbose=0;

	switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0)
            eprintf("usage: xparkeyval [-v]\n");
        verbose = 1;
        break;
    default:
        eprintf("usage: xparkeyval [-v]\n");
    }

	double x=1.0, y=2.0, z=3.0;
	ParKeyVal *pkv = NULL;

	pkv = ParKeyVal_add(pkv, "y", &y);
	pkv = ParKeyVal_add(pkv, "x", &x);
	pkv = ParKeyVal_add(pkv, "z", &z);

	double *ptr;
	ptr = ParKeyVal_get(pkv, "x");
	assert(ptr == &x);
	assert(*ptr == 1.0);

	ptr = ParKeyVal_get(pkv, "y");
	assert(ptr == &y);
	assert(*ptr == 2.0);
		   
	ptr = ParKeyVal_get(pkv, "z");
	assert(ptr == &z);
	assert(*ptr == 3.0);

	assert(NULL == ParKeyVal_get(pkv, "nonexistent"));

	if(verbose)
		ParKeyVal_print(pkv, stdout);

	ParKeyVal_free(pkv);

	unitTstResult("ParKeyVal", "OK");
    return 0;
}
