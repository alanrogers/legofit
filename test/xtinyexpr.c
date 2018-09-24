/**
 * @file xptrset.c
 * @author Alan R. Rogers
 * @brief Unit test for PtrSet
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "tinyexpr.h"
#include "misc.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>

#ifdef NDEBUG
#  error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {
    int         verbose = 0;

    if(argc==2 && 0==strcmp("-v", argv[1]))
        verbose = 1;
    else if(argc != 1) {
        fprintf(stderr,"usage: xtinyexpr [-v]\n");
        exit(1);
    }

    int status;
    double w=1.0, x=2.0, y=3.0, z=4.0;
    te_variable *pars = NULL;

    pars = te_variable_push(pars, "w", &w);
    pars = te_variable_push(pars, "x", &x);
    pars = te_variable_push(pars, "y", &y);
    pars = te_variable_push(pars, "z", &z);

    char formula[100];
    sprintf(formula, "%s", "w + x*y - z");

    te_expr *expr = te_compile(formula, pars, &status);
    printf("status=%d\n", status);
    if(expr == NULL) {
        fprintf(stderr,"%s:%d: parse error\n", __FILE__,__LINE__);
        fprintf(stderr,"  %s\n", formula);
        fprintf(stderr,"  %*s^\nError near here\n", status-1, "");
        exit(EXIT_FAILURE);
    }

    int len=100;
    const double *dep[len];
    len = te_dependencies(expr, len, dep);
    assert(len==4);
    int gotw=0, gotx=0, goty=0, gotz=0;
    for(int i=0; i<4; ++i) {
        if(dep[i] == &w)
            ++gotw;
        else if(dep[i] == &x)
            ++gotx;
        else if(dep[i] == &y)
            ++goty;
        else if(dep[i] == &z)
            ++gotz;
    }
    assert(gotw==1);
    assert(gotx==1);
    assert(goty==1);
    assert(gotz==1);

    double val = te_eval(expr);
    assert(val == w + x*y - z);

    if(verbose)
        te_print(expr, stdout);
    
    te_free(expr);

    unitTstResult("tinyexpr", "OK");
    return 0;
}
    
