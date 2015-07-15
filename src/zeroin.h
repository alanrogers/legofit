#ifndef LDPSIZ_ZEROIN
#define LDPSIZ_ZEROIN
int zeroin(double *root, double a, double b, double (*f) (double, void *),
           double tol, void *data, int verbose);
#endif
