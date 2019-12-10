/**
 * @file xmatcoal.c
 * @author Alan R. Rogers
 * @brief Test matcoal.c.
 * @copyright Copyright (c) 2019, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "matcoal.h"
#include <stdio.h>
#include <math.h>

int main(void) {

    int i, nsamples = 5;
    MatCoal_initExterns(nsamples);
    MatCoal_printAll(stdout);

    double ans[nsamples];
    double v = 0.01;
    int dim=2;

    printf("dim=%d\n", dim);
    MatCoal_project(dim, ans, v);
    printf("x[%lf]:", v);
    for(i=0; i<dim; ++i)
        printf(" %lg", ans[i]);
    putchar('\n');

    MatCoal_ciLen(dim, ans, v);
    printf("len[%lf]:", v);
    for(i=0; i<dim; ++i)
        printf(" %lg", ans[i]);
    fputs("\n\n", stdout);

    dim = nsamples-1;
    printf("dim=%d\n", dim);
    MatCoal_project(dim, ans, v);
    printf("x[%lf]:", v);
    for(i=0; i<dim; ++i)
        printf(" %lg", ans[i]);
    putchar('\n');

    MatCoal_ciLen(dim, ans, v);
    printf("len[%lf]:", v);
    for(i=0; i<dim; ++i)
        printf(" %lg", ans[i]);
    fputs("\n\n", stdout);

    v = 0.0;
    dim=2;

    printf("dim=%d\n", dim);
    MatCoal_project(dim, ans, v);
    printf("x[%lf]:", v);
    for(i=0; i<dim; ++i)
        printf(" %lg", ans[i]);
    putchar('\n');

    MatCoal_ciLen(dim, ans, v);
    printf("len[%lf]:", v);
    for(i=0; i<dim; ++i)
        printf(" %lg", ans[i]);
    fputs("\n\n", stdout);

    dim = nsamples-1;
    printf("dim=%d\n", dim);
    MatCoal_project(dim, ans, v);
    printf("x[%lf]:", v);
    for(i=0; i<dim; ++i)
        printf(" %lg", ans[i]);
    putchar('\n');

    MatCoal_ciLen(dim, ans, v);
    printf("len[%lf]:", v);
    for(i=0; i<dim; ++i)
        printf(" %lg", ans[i]);
    fputs("\n\n", stdout);

    v = INFINITY;
    dim=2;

    printf("dim=%d\n", dim);
    MatCoal_project(dim, ans, v);
    printf("x[%lf]:", v);
    for(i=0; i<dim; ++i)
        printf(" %lg", ans[i]);
    putchar('\n');

    MatCoal_ciLen(dim, ans, v);
    printf("len[%lf]:", v);
    for(i=0; i<dim; ++i)
        printf(" %lg", ans[i]);
    fputs("\n\n", stdout);

    dim = nsamples-1;
    printf("dim=%d\n", dim);
    MatCoal_project(dim, ans, v);
    printf("x[%lf]:", v);
    for(i=0; i<dim; ++i)
        printf(" %lg", ans[i]);
    putchar('\n');

    MatCoal_ciLen(dim, ans, v);
    printf("len[%lf]:", v);
    for(i=0; i<dim; ++i)
        printf(" %lg", ans[i]);
    fputs("\n\n", stdout);

    MatCoal_freeExterns();

    return 0;
}
