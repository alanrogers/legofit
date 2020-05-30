#ifndef ARR_MATCOAL_H
#define ARR_MATCOAL_H

#include <stdio.h>

void MatCoal_initExterns(long nsamp);
void MatCoal_freeExterns(void);
void MatCoal_eigenvals(int dim, double eig[dim], double v);
void MatCoal_project(int dim, double ans[dim], double eig[dim]);
void MatCoal_ciLen(int dim, double ans[dim], double eig[dim]);
void MatCoal_printAll(FILE *fp);

#endif
