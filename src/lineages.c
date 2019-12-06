#include "typedefs.h"
#include "misc.h"
#include "lblndx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct Lineages {
    struct Lineages *next;
    double prob;
    int dim; // allocated dimension of id
    tipId_t id[1];
} Lineages;

Lineages *Lineages_new(Lineages *next, int dim, tipId_t id[dim]);
Lineages *Lineages_coalesce(Lineages *self);
void Lineages_print(Lineages *self, FILE *fp);

Lineages *Lineages_new(Lineages *next, int dim, tipId_t id[dim]) {
    size_t size = sizeof(Lineages) + (dim-1)*sizeof(tipId_t);
    Lineages *new = malloc(size);
    CHECKMEM(new);

    new->next = next;
    new->prob = 0.0;
    new->dim = dim;
    memcpy(new->id, id, dim * sizeof(id[0]));
    return new;
}

// Return a list containing all possible states after
// a single coalescent event;
Lineages *Lineages_coalesce(Lineages *self) {
    Lineages *new = NULL;
    while(self) {
        for(int i=0; i < self->dim - 1; ++i) {
            for(int j=i+1; j < self->dim; ++j) {
                int k, dim2 = self->dim-1;
                tipId_t id2[dim2];
                for(k=0; k<i; ++k)
                    id2[k] = self->id[k];
                id2[i] = self->id[i] | self->id[j];
                for(k=i+1; k < j; ++k)
                    id2[k] = self->id[k];
                for(k=j; k < dim2; ++k)
                    id2[k] = self->id[k+1];
                qsort(id2, dim2, sizeof(id2[0]), compare_tipId);
                new = Lineages_new(new, dim2, id2);
                CHECKMEM(new);
            }
        }
        self = self->next;
    }
    return new;
}

void Lineages_print(Lineages *self, FILE *fp) {
    if(self==NULL)
        return;
    for(int i=0; i < self->dim; ++i)
        fprintf(fp, "%u ", self->id[i]);
    putc('\n', fp);
    Lineages_print(self->next, fp);
}

int main(void) {
    tipId_t id[] = {1, 2, 4, 8, 16};
    int dim = sizeof(id)/sizeof(id[0]);

    Lineages *lin0 = Lineages_new(NULL, dim, id);
    printf("0 coalescence\n");
    Lineages_print(lin0, stdout);

    Lineages *lin1 = Lineages_coalesce(lin0);
    printf("1 coalescence\n");
    Lineages_print(lin1, stdout);

    Lineages *lin2 = Lineages_coalesce(lin1);
    printf("2 coalescence\n");
    Lineages_print(lin2, stdout);

    Lineages *lin3 = Lineages_coalesce(lin2);
    printf("3 coalescence\n");
    Lineages_print(lin3, stdout);

    Lineages *lin4 = Lineages_coalesce(lin3);
    printf("4 coalescence\n");
    Lineages_print(lin4, stdout);

    
    return 0;
}
