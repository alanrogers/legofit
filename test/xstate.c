/**
 * @file xstate.c
 * @author Alan R. Rogers
 * @brief Test state.c.
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "typedefs.h"
#include "state.h"
#include "misc.h"
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>

#ifdef NDEBUG
#  error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {
	int verbose=0;

	switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0) {
            fprintf(stderr, "usage: xstate [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xstate [-v]\n");
        exit(EXIT_FAILURE);
    }

    const char *fname = "xstate.tmp";
    const int npts=3, npar=2;
    int i, status;
    double x[npts][npar] = {{1.0, 2.0}, {3.0,4.0}, {5.0, 6.0}};
    State *s = State_new(npts, npar);
    CHECKMEM(s);
    for(i=0; i < npts; ++i) {
        status = State_setVector(s, i, npar, x[i]);
        switch(status) {
        case 0:
            break;
        case EINVAL:
            fprintf(stderr,"%s:%d: Dimension mismatch in State_setVector\n",
                    __FILE__,__LINE__);
            exit(1);
        default:
            fprintf(stderr,"%s:%d: Unknown error\n", __FILE__,__LINE__);
            exit(1);
        }
    }

    FILE *fp = fopen(fname, "w");
    assert(fp);
    status = State_print(s, fp);
    switch(status) {
    case 0:
        break;
    case EIO:
        fprintf(stderr,"%s:%d: can't write to file\n", __FILE__,__LINE__);
        exit(1);
    default:
        fprintf(stderr,"%s:%d: Unknown error\n", __FILE__,__LINE__);
        exit(1);
    }

    State_free(s);

    fp = fopen(fname, "r");
    assert(fp);
    s = State_read(fp);
    CHECKMEM(s);
    fclose(fp);
    unlink(fname);

    double y[npar];
    for(i=0; i<npts; ++i) {
        State_getVector(s, i, npar, y);
        assert(0 == memcmp(y, x[i], npar*sizeof(y[0])));
    }

    unitTstResult("State", "OK");
    return 0;
}
