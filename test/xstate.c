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
    const char *fname2 = "xstate2.tmp";

    NameList *list=NULL;
    list = NameList_append(list, fname);
    list = NameList_append(list, fname2);
    assert(2 == NameList_size(list));

    const int npts=3, npar=2;
    int i, j, status;
    const char *parname[npar] = {"var1", "var2"};
    double x1[npts][npar] = {{1.0, 2.0}, {1.5, 3.5}, {2.0, 2.5}};
    double x2[npts][npar] = {{1.8, 2.2}, {1.7, 3.1}, {1.9, 2.8}};
    double c1[npts] = {0.01, 0.02, 0.03};
    double c2[npts] = {0.04, 0.05, 0.06};
    State *s = State_new(npts, npar);
    State *s2 = State_new(npts, npar);
    CHECKMEM(s);
    assert(npts == State_npoints(s));
    assert(npar == State_nparameters(s));
    assert(npts == State_npoints(s2));
    assert(npar == State_nparameters(s2));

    // set parameter names
    for(j=0; j < npar; ++j) {
        State_setName(s, j, parname[j]);
        State_setName(s2, j, parname[j]);
    }

    // set parameter values
    for(i=0; i < npts; ++i) {
        State_setVector(s, i, npar, x1[i]);
        State_setCost(s, i, c1[i]);
        State_setVector(s2, i, npar, x2[i]);
        State_setCost(s2, i, c2[i]);
    }

    // write old-type state file
    FILE *fp = fopen(fname, "w");
    assert(fp);
    fprintf(fp,"%d %d\n", npts, npar);
    for(i=0; i<npts; ++i) {
        fprintf(fp, "%0.18lf", c1[i]);
        for(j=0; j < npar; ++j)
            fprintf(fp, " %0.18lf", x1[i][j]);
        putc('\n', fp);
    }
    fclose(fp);

    // write new-type state file
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
    fclose(fp);

    State_free(s2);

    // read old-type file
    s = State_read(fname, npar, parname);
    CHECKMEM(s);

    // read new-type file
    s2 = State_read(fname2, npar, parname);
    CHECKMEM(s2);

    // check the two State objects
    double y[npar];
    for(i=0; i<npts; ++i) {
        State_getVector(s, i, npar, y);
        assert(0 == memcmp(y, x1[i], npar*sizeof(y[0])));
        assert(c1[i] == State_getCost(s, i));

        State_getVector(s2, i, npar, y);
        assert(0 == memcmp(y, x2[i], npar*sizeof(y[0])));
        assert(c2[i] == State_getCost(s2, i));
    }
    State_free(s);
    State_free(s2);

    // Read npts points spread across the two state files
    s = State_readList(list, npts, npar, parname);
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
