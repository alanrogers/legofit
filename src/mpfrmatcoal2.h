#ifndef ARR_MPFRMATCOAL2_H
#define ARR_MPFRMATCOAL2_H

#include <stdio.h>

void MpfrMatCoal_initExterns(long nsamp);
void MpfrMatCoal_freeExterns(void);
void MpfrMatCoal_project(int dim, double ans[dim], double v);
void MpfrMatCoal_ciLen(int dim, double ans[dim], double v);
void MpfrMatCoal_printAll(FILE *fp);

#endif
