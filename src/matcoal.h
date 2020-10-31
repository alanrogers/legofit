#ifndef ARR_MATCOAL_H
#define ARR_MATCOAL_H

#include <stdio.h>

int  MatCoal_nSamples(void);
void MatCoal_initExterns(long nsamp);
void MatCoal_freeExterns(void);
void MatCoal_eigenvals(int dim, long double eig[dim], long double v);
void MatCoal_project(int dim, long double ans[dim], long double eig[dim]);
void MatCoal_ciLen(int dim, long double ans[dim], long double eig[dim]);
void MatCoal_printAll(FILE *fp);

#endif
