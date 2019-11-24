#include "typedefs.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef struct {
    double prob;
    int dim; // allocated dimension of id
    int n;   // number of ids stored
    tipId_t id[1];
} Lineages;

Lineages *Lineages_new(int dim, tipId_t id[dim]);

Lineages *Lineages_new(int dim, tipId_t id[dim]) {
    size_t size = sizeof(Lineages) + (dim-1)*sizeof(tipId_t);
    Lineages *new = malloc(size);
    CHECKMEM(new);

    new->prob = 0.0;
    new->n = 0;
    new->dim = dim;
    memcpy(new->id, id, dim * sizeof(dim[0]));
    return new;
}

int main(void) {
    unsigned id0[5] = {1, 2, 4, 8, 16};
    
}
