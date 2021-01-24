#ifndef ARR_MPFRMATCOAL_H
#define ARR_MPFRMATCOAL_H

#include <stdio.h>

int  MpfrMatCoal_isInitialized(void);
void MpfrMatCoal_initExterns(long nsamp);
void MpfrMatCoal_freeExterns(void);
void MpfrMatCoal_project(int dim, long double ans[dim], long double v);
void MpfrMatCoal_ciLen(int dim, long double ans[dim], long double v);
void MpfrMatCoal_printAll(FILE *fp);

#endif
