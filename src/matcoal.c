#include "rational.h"
#include <stdio.h>
#include <string.h>

// Number of haploid samples in the data. 0 until initialized
static int nsamples = 0;

// Array of pointers to matrices of colum
static double **cvec = NULL;
