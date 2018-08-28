/**
 * @file strparmap.c
 * @brief Map parameter names to parameter structures.
 * 
 * This class uses a left-leaning binary tree to build a map
 * that associates parameter names with structures that include
 * a pointer to the current value, the lower and upper bounds, and
 * the type of the parameter. The interface allows parameters to be
 * inserted and later looked up by name.
 *
 * Translated from Java code at
 * http://www.cs.princeton.edu/~rs/talks/LLRB/LLRB.pdf 
 */
#include "param.h"
#include "strparmap.h"
#include "misc.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#define RED 1
#define BLACK 0
#define ISRED(node) ((node)!=NULL && (node)->color==RED)

struct StrParMap {
    Param *par; // locally owned
    int color;
    StrParMap *left, *right;
};

StrParMap *StrParMap_new(Param *par);
StrParMap *StrParMap_insert_r(StrParMap *h, Param *par);
StrParMap *StrParMap_rotateLeft(StrParMap *h);
StrParMap *StrParMap_rotateRight(StrParMap *h);
void       StrParMap_flipColors(StrParMap *h);

void StrParMap_free(StrParMap *h) {
    if(h == NULL)
        return;
    StrParMap_free(h->left);
    StrParMap_free(h->right);
    Param_free(h->par);
    free(h);
}

const char *StrParMap_key(StrParMap *h) {
    return h->par->name;
}

StrParMap *StrParMap_new(Param *par) {
    StrParMap *self = malloc(sizeof(StrParMap));
    if(self==NULL) {
        fprintf(stderr,"%s:%d: can't allocate memory\n",
                __FILE__, __LINE__);
        exit(1);
    }
    self->par = par;
    self->color = RED;
    return self;
}

Param *StrParMap_search(StrParMap *root, const char *key) {
    StrParMap *x = root;
    while(x != NULL) {
        int diff = strcmp(key, x->par->name);
        if(diff == 0)
            return x->par;
        else if(diff < 0)
            x = x->left;
        else
            x = x->right;
    }
    return NULL;
}

StrParMap *StrParMap_insert(StrParMap *root, Param *par) {
    root = StrParMap_insert_r(root, par);
    root->color = BLACK;
    return root;
}

StrParMap *StrParMap_insert_r(StrParMap *h, Param *par) {
    if(h == NULL)
        return StrParMap_new(par);
    if(ISRED(h->left) && ISRED(h->right))
        StrParMap_flipColors(h);
    int diff = strcmp(par->name, h->par->name);
    if(diff == 0) {
        fprintf(stderr,"%s:%d: duplicate parameter name: %s\n",
                __FILE__,__LINE__, par->name);
        exit(EXIT_FAILURE);
    }
    if( diff < 0 )
        h->left = StrParMap_insert_r(h->left, par);
    else {
        assert(diff > 0);
        h->right = StrParMap_insert_r(h->right, par);
    }
    if(ISRED(h->right) && !ISRED(h->left))
        h = StrParMap_rotateLeft(h);
    if(ISRED(h->left) && ISRED(h->left->left))
        h = StrParMap_rotateRight(h);
    return h;
}

StrParMap *StrParMap_rotateLeft(StrParMap *h) {
    StrParMap *x = h->right;
    h->right = x->left;
    x->left = h;
    x->color = h->color;
    h->color = RED;
    return x;
}

StrParMap *StrParMap_rotateRight(StrParMap *h) {
    StrParMap *x = h->left;
    h->left = x->right;
    x->right = h;
    x->color = h->color;
    h->color = RED;
    return x;
}

void StrParMap_flipColors(StrParMap *h) {
    h->color = !h->color;
    h->left->color = !h->left->color;
    h->right->color = !h->right->color;
}

#ifdef TEST

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

#include <string.h>

void StrParMap_print(StrParMap *h, FILE *fp, int indent);

void StrParMap_print(StrParMap *h, FILE *fp, int indent) {
    if(h == NULL)
        return;
    StrParMap_print(h->left, fp, indent+1);
    int i;
    for(i=0; i<indent; ++i)
        fprintf(fp,"%d", i);
    fprintf(fp, "[%s, %p, %s]\n",
            h->par->name,
            h->par->valptr,
            h->color==RED ? "red" : "black");
    StrParMap_print(h->right, fp, indent+1);
}

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
#endif
