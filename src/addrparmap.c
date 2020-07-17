/**
 * @file addrparmap.c
 * @brief Map parameter addresses to parameter structures.
 * 
 * This class uses a left-leaning binary tree to build a map
 * that associates parameter addresses with structures that include
 * a pointer to the current value, the lower and upper bounds, and
 * the type of the parameter. The interface allows parameters to be
 * inserted and later looked up by address.
 *
 * Translated from Java code at
 * http://www.cs.princeton.edu/~rs/talks/LLRB/LLRB.pdf 
 */
#include "addrparmap.h"
#include "param.h"
#include "misc.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#define RED 1
#define BLACK 0
#define ISRED(node) ((node)!=NULL && (node)->color==RED)

struct AddrParMap {
    Param *par; // not locally owned
    int color;
    AddrParMap *left, *right;
};

AddrParMap *AddrParMap_new(Param *par);
AddrParMap *AddrParMap_insert_r(AddrParMap *h, Param *par);
AddrParMap *AddrParMap_rotateLeft(AddrParMap *h);
AddrParMap *AddrParMap_rotateRight(AddrParMap *h);
void AddrParMap_flipColors(AddrParMap *h);

void AddrParMap_free(AddrParMap *h) {
    if(h == NULL)
        return;
    AddrParMap_free(h->left);
    AddrParMap_free(h->right);
    free(h);
}

AddrParMap *AddrParMap_new(Param *par) {
    AddrParMap *self = malloc(sizeof(AddrParMap));
    CHECKMEM(self);
    self->par = par;  // no copy because par is not locally owned
    self->color = RED;
    self->left = self->right = NULL;
    return self;
}

Param *AddrParMap_search(AddrParMap *root, const double *valptr) {
    AddrParMap *x = root;
    while(x != NULL) {
        if(valptr == &x->par->value)
            return x->par;
        else if(valptr < &x->par->value)
            x = x->left;
        else
            x = x->right;
    }
    return NULL;
}

AddrParMap *AddrParMap_insert(AddrParMap *root, Param *par) {
    root = AddrParMap_insert_r(root, par);
    root->color = BLACK;
    return root;
}

AddrParMap *AddrParMap_insert_r(AddrParMap *h, Param *par) {
    if(h == NULL)
        return AddrParMap_new(par);
    if(ISRED(h->left) && ISRED(h->right))
        AddrParMap_flipColors(h);
    if(&par->value == &h->par->value) {
        fprintf(stderr,"%s:%d: duplicate value pointer: %s\n",
                __FILE__,__LINE__, par->name);
        exit(EXIT_FAILURE);
    }
    if( &par->value < &h->par->value )
        h->left = AddrParMap_insert_r(h->left, par);
    else {
        assert(&par->value > &h->par->value);
        h->right = AddrParMap_insert_r(h->right, par);
    }
    if(ISRED(h->right) && !ISRED(h->left))
        h = AddrParMap_rotateLeft(h);
    if(h->left && ISRED(h->left) && ISRED(h->left->left))
        h = AddrParMap_rotateRight(h);
    return h;
}

AddrParMap *AddrParMap_rotateLeft(AddrParMap *h) {
    AddrParMap *x = h->right;
    h->right = x->left;
    x->left = h;
    x->color = h->color;
    h->color = RED;
    return x;
}

AddrParMap *AddrParMap_rotateRight(AddrParMap *h) {
    AddrParMap *x = h->left;
    h->left = x->right;
    x->right = h;
    x->color = h->color;
    h->color = RED;
    return x;
}

void AddrParMap_flipColors(AddrParMap *h) {
    h->color = !h->color;
    h->left->color = !h->left->color;
    h->right->color = !h->right->color;
}

#ifdef TEST

#  ifdef NDEBUG
#    error "Unit tests must be compiled without -DNDEBUG flag"
#  endif

#include "network.h"
#include <string.h>

void AddrParMap_print(AddrParMap *h, FILE *fp, int indent);

void AddrParMap_print(AddrParMap *h, FILE *fp, int indent) {
    if(h == NULL)
        return;
    AddrParMap_print(h->left, fp, indent+1);
    int i;
    for(i=0; i<indent; ++i)
        fprintf(fp,"%d", i);

    // Printing address of value, because that's the sort key.
    fprintf(fp, "[%p, %s, %s]\n",
            &h->par->value,
            h->par->name,
            h->color==RED ? "red" : "black");
    AddrParMap_print(h->right, fp, indent+1);
}

int main(int argc, char **argv) {

    int         verbose = 0;

    if(argc > 1) {
        if(argc != 2 || 0 != strcmp(argv[1], "-v")) {
            fprintf(stderr, "usage: xaddrparmap [-v]\n");
            exit(EXIT_FAILURE);
        }
        verbose = 1;
    }

    Network_init(SIM);

    double a[] = {0.0, 1.0, 2.0, 0.5};

    Param *parvec[4];
    parvec[0] = Param_new("par0", a[0], 0.0, 1.0, FREE, NULL);
    parvec[1] = Param_new("par1", a[1], -INFINITY, INFINITY,
                          CONSTRAINED, "3 + 4*y");
    parvec[2] = Param_new("par2", a[2], 1.0, 1e6, FIXED, NULL);
    parvec[3] = Param_new("par3", a[3], -1.0, 1.0, CONSTRAINED,
                          "y - z");

    AddrParMap *root = NULL;
    root = AddrParMap_insert(root, parvec[0]);
    root = AddrParMap_insert(root, parvec[1]);
    root = AddrParMap_insert(root, parvec[2]);
    root = AddrParMap_insert(root, parvec[3]);

    if(verbose)
        AddrParMap_print(root, stdout, 0);

    // Search for nodes that exist. Each search should succeed.
    Param *par;

    par = AddrParMap_search(root, &(parvec[0]->value));
    assert(par == parvec[0]);
    assert(strcmp(par->name, "par0") == 0);

    par = AddrParMap_search(root, &(parvec[1]->value));
    assert(par == parvec[1]);
    assert(strcmp(par->name, "par1") == 0);
    
    par = AddrParMap_search(root, &(parvec[2]->value));
    assert(par == parvec[2]);
    assert(strcmp(par->name, "par2") == 0);

    par = AddrParMap_search(root, &(parvec[3]->value));
    assert(par == parvec[3]);
    assert(strcmp(par->name, "par3") == 0);

    // Search for nodes that don't exist. Each search should fail.
    par = AddrParMap_search(root, 1lu + &(parvec[3]->value));
    assert(par == NULL);

    AddrParMap_free(root);
    for(int i=0; i < 4; ++i)
        free(parvec[i]);
    printf("%-26s %s\n", "AddrParMap", "OK");
}
#endif
