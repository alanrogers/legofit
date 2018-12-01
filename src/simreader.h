#ifndef SIMREADER_H
#define SIMREADER_H

#include "typedefs.h"
#include <stdio.h>

// constructor
SimReader *SimReader_new(FILE *fp);

// destructor
void SimReader_free(SimReader *self);

// Rewind input and reset chr and nucpos. Doesn't work
// if input is stdin.
int SimReader_rewind(SimReader *self);

// Move SimReader to next nucleotide site.
int SimReader_next(SimReader *self);

// Return current chromosome.
unsigned SimReader_chr(SimReader *self);

// Return the dimension of the array of samples
int SimReader_sampleDim(SimReader *self);

// Return number of samples from population i.
int SimReader_nsamples(SimReader *self, int i);

// Return frequency of derived allele in sample from population i.
double SimReader_daf(SimReader *self, int i);

// Return pointer to label of population i
const char *SimReader_lbl(SimReader *self, int i);
#endif
