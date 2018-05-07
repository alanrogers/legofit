#include "wraparound.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
int main(int argc, char **argv) {
    int verbose = 0;
    unsigned i, u;

	switch (argc) {
    case 1:
        break;
    case 2:
        if(strncmp(argv[1], "-v", 2) != 0) {
            fprintf(stderr, "usage: xwraparound [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xwraparound [-v]\n");
        exit(EXIT_FAILURE);
    }

    unsigned totsize = 8;
    Wraparound *w = Wraparound_new(totsize);

    assert(0 == Wraparound_size(w));
    Wraparound_push(w, 0u);
    assert(0 != Wraparound_size(w));
    assert(0 == Wraparound_pop(w));
    assert(0 == Wraparound_size(w));
    for(u=0u; u < 2*totsize; ++u)
        Wraparound_push(w, u);
    assert(totsize == Wraparound_size(w));

    for(i=0; i < totsize/2; ++i)
        (void) Wraparound_pop(w);
    assert(0 != Wraparound_size(w));
    assert(totsize != Wraparound_size(w));
    
    for(u=0u; u < totsize/2u; ++u)
        Wraparound_push(w, u);
    assert(totsize == Wraparound_size(w));

    if(verbose) {
        printf("Popping Wraparound:\n");
        Wraparound_print(w, stdout);
        while(0 != Wraparound_size(w))
            printf("%u\n", Wraparound_pop(w));
    }

    Wraparound_free(w);
    
    printf("Wraparound OK\n");
    return 0;
}
