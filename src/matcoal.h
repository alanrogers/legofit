#ifndef ARR_MATCOAL_H
#define ARR_MATCOAL_H

int MatCoal_initExterns(int n);
void MatCoal_freeExterns(void);
void MatCoal_project(int dim, double ans[dim], double v);
void MatCoal_ciLen(int dim, double ans[dim], double v);

#endif
