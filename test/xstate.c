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

    StateFile_t statefiletype = OLD;

    const char *fname = "xstate.tmp";
    const char *fname2 = "xstate2.tmp";

    NameList *list=NULL;
    list = NameList_append(list, fname);
    list = NameList_append(list, fname2);
    assert(2 == NameList_size(list));

    const int npts=3, npar=2;
    int i, status;
    double x1[npts][npar] = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
    double x2[npts][npar] = {{7.0, 8.0}, {9.0, 10.0}, {11.0, 12.0}};
    double c1[npts] = {0.01, 0.02, 0.03};
    double c2[npts] = {0.04, 0.05, 0.06};
    State *s = State_new(npts, npar);
    State *s2 = State_new(npts, npar);
    CHECKMEM(s);
    assert(npts == State_npoints(s));
    assert(npar == State_nparameters(s));
    for(i=0; i < npts; ++i) {
        State_setVector(s, i, npar, x1[i]);
        State_setCost(s, i, c1[i]);
        State_setVector(s2, i, npar, x2[i]);
        State_setCost(s2, i, c2[i]);
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
        fprintf(stderr,"%s:%d: Unknown error %d\n", __FILE__,__LINE__,
                status);
        exit(1);
    }
    fclose(fp);

    fp = fopen(fname2, "w");
    assert(fp);
    status = State_print(s2, fp);
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

    State_free(s2);

    fclose(fp);
    fp = fopen(fname, "r");
    assert(fp);
    s = State_read(fp, statefiletype);
    CHECKMEM(s);
    fclose(fp);

    double y[npar];
    for(i=0; i<npts; ++i) {
        State_getVector(s, i, npar, y);
        assert(0 == memcmp(y, x1[i], npar*sizeof(y[0])));
        assert(c1[i] == State_getCost(s, i));
    }
    State_free(s);

    s = State_readList(list, npts, gptree, statefiletype);
    if(verbose)
        State_print(s, stderr);
    for(i=0; i<npts; ++i) {
        State_getVector(s, i, npar, y);
        if(i < 2) {
            assert(0 == memcmp(y, x1[i], npar*sizeof(y[0])));
            assert(c1[i] == State_getCost(s, i));
        }else{
            assert(0 == memcmp(y, x2[i-2], npar*sizeof(y[0])));
            assert(c2[i-2] == State_getCost(s, i));
        }
    }
    State_free(s);

    NameList_free(list);
    unitTstResult("NameList", "OK");
    unitTstResult("State", "OK");

    unlink(fname);
    unlink(fname2);
    return 0;
}
