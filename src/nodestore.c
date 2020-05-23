/**
 * @file nodestore.c
 * @author Alan R. Rogers
 * @brief Allocate nodes from an array
 *
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "misc.h"
#include "nodestore.h"

/// This structure allows you to allocate PopNode objects in an array
/// and then dole them out one at a time via calls to NodeStore_alloc.
/// This keeps all the PopNode objects together in memory and may
/// reduce page faults.
struct NodeStore {
    size_t curr, end, elsize;
    void *v;              // not locally owned
};

/// Allocate a new NodeStore, which provides an interface
/// for getting PopNode objects, one at a time, out
/// of a previously-allocated array v.
/// @param[in] len number of PopNode objects in array.
/// @param[in] v array of PopNode objects
/// @return newly-allocated NodeStore.
NodeStore  *NodeStore_new(unsigned len, size_t elsize, void * v) {
    NodeStore  *self = malloc(sizeof(NodeStore));
    CHECKMEM(self);

    self->v = v;
    self->curr = (size_t) self->v;
    self->elsize = elsize;
    self->end = self->curr + len * self->elsize;

    return self;
}

/// Destructor for NodeStore
void NodeStore_free(NodeStore * self) {
    // Does not free self->v
    free(self);
}

/// Return a pointer to an unused element within NodeStore. Return
/// NULL if none are left.
void *NodeStore_alloc(NodeStore * self) {
    if(self->curr >= self->end) {
        assert(self->curr == self->end);
        return NULL;
    }
    void *new = (void *) self->curr;
    self->curr += self->elsize;
    return new;
}

#ifdef TEST

#  include <string.h>
#  include <assert.h>
#  include <time.h>

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

int main(int argc, char **argv) {

    int         verbose = 0;

    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xnodestore [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    int         nseg = 10;
    PopNode     v[nseg];
    NodeStore  *ns = NodeStore_new(nseg, sizeof(v[0]), v);
    CHECKMEM(ns);

    PopNode    *node;
    for(int i=0; i < nseg; ++i) {
        node = NodeStore_alloc(ns);
        assert(node = v+i);
    }
    node = NodeStore_alloc(ns);
    assert(node == NULL);

    NodeStore_free(ns);
    unitTstResult("NodeStore", "OK");

    return 0;
}
#endif
