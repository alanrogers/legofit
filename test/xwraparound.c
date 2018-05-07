#include "ringbuf.h"
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
            fprintf(stderr, "usage: xringbuf [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
        break;
    default:
        fprintf(stderr, "usage: xringbuf [-v]\n");
        exit(EXIT_FAILURE);
    }

    RingBuf *rb = RingBuf_new();

    assert(RingBuf_empty(rb));
    assert(!RingBuf_full(rb));
    RingBuf_pushHead(rb, 0u);
    assert(!RingBuf_empty(rb));
    assert(!RingBuf_full(rb));
    assert(0 == RingBuf_popTail(rb));
    assert(RingBuf_empty(rb));
    assert(!RingBuf_full(rb));
    for(u=0u; u < 2*RBUF_SIZE; ++u)
        RingBuf_pushHead(rb, u);
    assert(RingBuf_full(rb));

    for(i=0; i < RBUF_SIZE/2; ++i)
        (void) RingBuf_popTail(rb);
    assert(!RingBuf_full(rb));
    
    for(u=0u; u < RBUF_SIZE/2u; ++u)
        RingBuf_pushHead(rb, u);
    assert(RingBuf_full(rb));

    if(verbose) {
        printf("Popping RingBuf:\n");
        RingBuf_print(rb, stdout);
        while(!RingBuf_empty(rb)) {
            printf("%u\n", RingBuf_popTail(rb));
        }
    }

    RingBuf_free(rb);
    
    printf("RingBuf OK\n");
    return 0;
}
