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
#include "gptree.h"
#include "parstore.h"
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>

#ifdef NDEBUG
#  error "Unit tests must be compiled without -DNDEBUG flag"
#endif

const char *tstInput =
    "# this is a comment\n"
    "time fixed  T0=0\n"
    "time free   Tc=1\n"
    "time free   Tab=3\n"
    "time fixed  Tabc=5.5\n"
    "twoN fixed  twoNa=100\n"
    "twoN fixed  twoNb=123\n"
    "twoN fixed  twoNc=213.4\n"
    "twoN fixed  twoNbb=32.1\n"
    "twoN fixed  twoNab=222\n"
    "twoN fixed  twoNabc=1.2e2\n"
    "mixFrac fixed Mc=0.02\n"
    "segment a   t=T0     twoN=twoNa    samples=1\n"
    "segment b   t=T0     twoN=twoNb    samples=1\n"
    "segment c   t=Tc     twoN=twoNc    samples=1\n"
    "segment bb  t=Tc     twoN=twoNbb\n"
    "segment ab  t=Tab    twoN=twoNab\n"
    "segment abc t=Tabc   twoN=twoNabc\n"
    "mix    b  from bb + Mc * c\n"
    "derive a  from ab\n"
    "derive bb from ab\n" "derive ab from abc\n" "derive c  from abc\n";

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

    const char *lgoname = "xstate.lgo";
    const char *fname = "xstate.tmp";
    const char *fname2 = "xstate2.tmp";
    FILE       *fp = fopen(lgoname, "w");
    fputs(tstInput, fp);
    fclose(fp);

    Bounds      bnd = {
        .lo_twoN = 0.0,
        .hi_twoN = 1e7,
        .lo_t = 0.0,
        .hi_t = HUGE_VAL
    };

    GPTree     *gptree = GPTree_new(lgoname, bnd);
    assert(gptree);

    NameList *list=NULL;
    list = NameList_append(list, fname);
    list = NameList_append(list, fname2);
    assert(2 == NameList_size(list));

    const int npts=3, npar=2;
    int i, j, status;
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
        State_setName(s, j, GPTree_getNameFree(gptree, j));
        State_setName(s2, j, GPTree_getNameFree(gptree, j));
    }

    // set parameter values
    for(i=0; i < npts; ++i) {
        State_setVector(s, i, npar, x1[i]);
        State_setCost(s, i, c1[i]);
        State_setVector(s2, i, npar, x2[i]);
        State_setCost(s2, i, c2[i]);
    }

    // write old-type state file
    fp = fopen(fname, "w");
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
    fp = fopen(fname, "r");
    assert(fp);
    s = State_read(fp);
    CHECKMEM(s);
    fclose(fp);

    // read new-type file
    fp = fopen(fname2, "r");
    assert(fp);
    s2 = State_read(fp);
    CHECKMEM(s2);
    fclose(fp);

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
    s = State_readList(list, npts, gptree);
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
    unlink(lgoname);
    return 0;
}
