#ifndef ARR_POINT_H
#define ARR_POINT_H

#include "typedefs.h"

void Point_setNPar(int npar);
Point *Point_new(void);
void Point_free(Point *self);
void Point_set(Point *self, double cost, double *par);
double Point_get(Point *self, double *par);

#endif
