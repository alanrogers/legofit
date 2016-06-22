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
	double xx=1.0, yy=2.0, zz=3.0;
	ParKeyVal *pkv = NULL;
	bool isfree;

	pkv = ParKeyVal_add(pkv, "y", &y, true);
	pkv = ParKeyVal_add(pkv, "x", &x, false);
	pkv = ParKeyVal_add(pkv, "z", &z, true);

	double *ptr;
	ptr = ParKeyVal_get(pkv, &isfree, "x");
	assert(ptr == &x);
	assert(isfree == false);
	assert(*ptr == 1.0);

	ptr = ParKeyVal_get(pkv, &isfree, "y");
	assert(ptr == &y);
	assert(isfree == true);
	assert(*ptr == 2.0);
		   
	ptr = ParKeyVal_get(pkv, &isfree, "z");
	assert(ptr == &z);
	assert(isfree == true);
	assert(*ptr == 3.0);

	assert(NULL == ParKeyVal_get(pkv, &isfree, "nonexistent"));
	assert(isfree == true);

	if(verbose)
		ParKeyVal_print(pkv, stdout);

    ParKeyVal *pkv2 = NULL;

	pkv2 = ParKeyVal_add(pkv2, "x", &xx, false);
	pkv2 = ParKeyVal_add(pkv2, "z", &zz, true);
	pkv2 = ParKeyVal_add(pkv2, "y", &yy, true);

    assert(ParKeyVal_equals(pkv, pkv2));

	ParKeyVal_free(pkv);
	ParKeyVal_free(pkv2);

	unitTstResult("ParKeyVal", "OK");
    return 0;
}
