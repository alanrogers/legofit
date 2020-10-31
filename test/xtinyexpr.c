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
#include "strptrmap.h"
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
    StrPtrMap *pars = StrPtrMap_new();

    status = StrPtrMap_insert(pars, "w", te_variable_new(&w));
    if(status) {
        fprintf(stderr, "%s:%d: duplicate insertion\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    status = StrPtrMap_insert(pars, "x", te_variable_new(&x));
    if(status) {
        fprintf(stderr, "%s:%d: duplicate insertion\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    
    status = StrPtrMap_insert(pars, "y", te_variable_new(&y));
    if(status) {
        fprintf(stderr, "%s:%d: duplicate insertion\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    status = StrPtrMap_insert(pars, "z", te_variable_new(&z));
    if(status) {
        fprintf(stderr, "%s:%d: duplicate insertion\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    char formula[100];
    sprintf(formula, "%s", "w + x*y - z");

    status = 0;
    te_expr *expr = te_compile(formula, pars, &status);
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

    sprintf(formula, "%s", "exp(x+log(y+z))");
    expr = te_compile(formula, pars, &status);
    if(expr == NULL) {
        fprintf(stderr,"%s:%d: parse error\n", __FILE__,__LINE__);
        fprintf(stderr,"  %s\n", formula);
        fprintf(stderr,"  %*s^\nError near here\n", status-1, "");
        exit(EXIT_FAILURE);
    }
    val = te_eval(expr);
    assert(val == exp(x+log(y+z)));
    
    te_free(expr);
    te_free_variables(pars);

    unitTstResult("tinyexpr", "OK");
    return 0;
}
    
