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

static void GET_PLOT_DATA(NODETYPE *self, FILE * fp,
                                         PtrQueue *time_zero,
                                         PtrQueue *dual_parents,
                                         PtrQueue *edge);

void PLOT(NODETYPE *root, SampNdx *sndx, FILE *fp) {
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

    if(PtrQueue_size(time_zero) > 1) {
        fputs("  { rank = same;", fp);
        while( (s=PtrQueue_pop(time_zero)) != NULL) {
            fprintf(fp, " %s", s);
            free(s);
        }
        fputs("};\n", fp);
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
