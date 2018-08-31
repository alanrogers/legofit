#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

#include "strparmap2.h"
#include "param.h"
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv) {

    int         verbose = 0;

    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xstrparmap [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    double a[] = {0.0, 1.0, 2.0, 3.0};

    StrParMap *root = NULL;
    root = StrParMap_insert(root, Param_new("par0", a+0, 0.0, 1.0, Free) );
    root = StrParMap_insert(root, Param_new("par1", a+1, -INFINITY, INFINITY,
                                       Constrained) );
    root = StrParMap_insert(root, Param_new("par2", a+2, 1.0, 1e6, Fixed) );
    root = StrParMap_insert(root, Param_new("par3", a+3, -1.0, 1.0, Gaussian) );

    if(verbose)
        StrParMap_print(root, stdout, 0);

    // Search for nodes that exist. Each search should succeed.
    Param *par;

    par = StrParMap_search(root, "par0");
    assert(par);
    assert(strcmp(par->name, "par0") == 0);

    par = StrParMap_search(root, "par1");
    assert(par);
    assert(strcmp(par->name, "par1") == 0);
    
    par = StrParMap_search(root, "par2");
    assert(par);
    assert(strcmp(par->name, "par2") == 0);

    par = StrParMap_search(root, "par3");
    assert(par);
    assert(strcmp(par->name, "par3") == 0);

    // Search for nodes that don't exist. Each search should fail.
    par = StrParMap_search(root, "par10");
    assert(par == NULL);

    StrParMap_free(root);
    printf("%-26s %s\n", "StrParMap", "OK");
}
