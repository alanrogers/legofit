/**
Generic code for functions that read a network to generate a file in
.dot format, which contains commands that produce a network graph. To
create a concrete version, include this into a .c file that defines
NODETYPE.
**/

// Macros for creating function names
#define INNER(CLASS, NAME) CLASS ## _ ## NAME
#define FUNC(CLASS, NAME) INNER(CLASS, NAME)

// Macros below expand to function names of form Foo_plot,
// Foo_get_plot_data, and Foo_unvisit, where "Foo" is the value of the
// NODETYPE macro.
#define PLOT FUNC(NODETYPE, plot)
#define GET_PLOT_DATA FUNC(NODETYPE, get_plot_data)
#define UNVISIT FUNC(NODETYPE, unvisit)

//#include "ptrqueue.h"

#include "sampndx.h"
#include "lblndx.h"

static void GET_PLOT_DATA(NODETYPE *self, FILE * fp,
                                         PtrQueue *time_zero,
                                         PtrQueue *dual_parents,
                                         PtrQueue *edge);

void PLOT(NODETYPE *root, SampNdx *sndx, LblNdx *lblndx, FILE *fp) {
    PtrQueue *time_zero = PtrQueue_new();
    PtrQueue *dual_parents = PtrQueue_new();
    PtrQueue *edge = PtrQueue_new();

    fprintf(fp, "digraph {\n");
    fprintf(fp, "  ratio = 1.5;\n");

    UNVISIT(root);
    GET_PLOT_DATA(root, fp, time_zero, dual_parents, edge);

    char *s;

    for(int i=0; i < SampNdx_size(sndx); ++i) {
        NODETYPE *node = SampNdx_get(sndx, i);
        fprintf(fp, "  %s [shape=square];\n", node->label);
    }

    // All nodes that occur at time zero should be on a single
    // horizontal rank. They should also appear left-to-right
    // in the order they occur in the .lgo file. To accomplish this,
    // the following paragraph generates a line of .dot input of the
    // following form:
    //
    //     { rank = same; x -> y -> z [style = invis]}
    //
    // where x, y, and z are the nodes at time 0 and they are listed
    // in the order of their appearance in the .lgo file. The string
    // "[style = invis]" instructs Graphviz not to print the arrows
    // specified here. Their only purpose is to specify the horizontal
    // order of x, y, and z.
    int nsamples = PtrQueue_size(time_zero);
    if(nsamples > 1) {
        fputs("  { rank = same;", fp);

        char *lbl[nsamples];
        tipId_t pat[nsamples];
        unsigned ord[nsamples];

        for(int i=0; i < nsamples; ++i) {
            lbl[i]=PtrQueue_pop(time_zero);
            pat[i] = LblNdx_getTipId(lblndx, lbl[i]);
        }

        orderpat(nsamples, ord, pat);

        for(int i = 0; i < nsamples-1; ++i) {
            fprintf(fp, " %s ->", lbl[ord[i]]);
        }
        fprintf(fp, " %s [style = invis]};\n", lbl[ord[nsamples-1]]);

        for(int i=0; i<nsamples; ++i)
            free(lbl[i]);
    }

    while( (s=PtrQueue_pop(dual_parents)) != NULL) {
        fprintf(fp, "  { rank = same; %s};\n", s);
        free(s);
    }
    
    while( (s=PtrQueue_pop(edge)) != NULL) {
        fprintf(fp, "  %s;\n", s);
        free(s);
    }

    fprintf(fp, "}\n");

    PtrQueue_free(time_zero);
    PtrQueue_free(dual_parents);
    PtrQueue_free(edge);
}

static void GET_PLOT_DATA(NODETYPE *self, FILE * fp,
                          PtrQueue *time_zero,
                          PtrQueue *dual_parents,
                          PtrQueue *edge) {

    char s[200] = {0};
    
    if(self->visited)
        return;
    self->visited = 1;

    if(self->start == 0) {
        PtrQueue_push(time_zero, strdup(self->label));
    }
    if(self->nparents == 2) {
        sprintf(s, "%s %s", self->parent[0]->label,
                self->parent[1]->label);
        PtrQueue_push(dual_parents, strdup(s));
    }
    if(self->nchildren == 1) {
        sprintf(s, "%s -> %s", self->label, self->child[0]->label);
        PtrQueue_push(edge, strdup(s));
    }else if(self->nchildren == 2) {
        sprintf(s, "%s -> {%s %s}", self->label,
                self->child[0]->label,
                self->child[1]->label);
        PtrQueue_push(edge, strdup(s));
    }

    for(int i = 0; i < self->nchildren; ++i)
        GET_PLOT_DATA(self->child[i], fp, time_zero, dual_parents, edge);
}

#undef INNER
#undef FUNC
#undef PLOT
#undef GET_PLOT_DATA
#undef UNVISIT
