/**
 * @file xuintqueue.c
 * @author Alan R. Rogers
 * @brief Test uintqueue.c.
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "uintqueue.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

int main(int argc, char **argv) {
    int         verbose = 0;
    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xuintqueue [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    UINTqueue *queue = NULL;
    queue = UINTqueue_push(queue, 1u);
    queue = UINTqueue_push(queue, 2u);
    assert(2 == UINTqueue_length(queue));
    unsigned x=0;
    queue = UINTqueue_pop(queue, &x);
    assert(1u == x);
    assert(1 == UINTqueue_length(queue));
    queue = UINTqueue_pop(queue, &x);
    assert(2u == x);
    assert(0 == UINTqueue_length(queue));

    queue = NULL;
    queue = UINTqueue_push(queue, 1u);
    queue = UINTqueue_push(queue, 2u);
    queue = UINTqueue_free(queue);
    assert(queue == NULL);
    unitTstResult("UINTqueue", "OK");

    return 0;
}
