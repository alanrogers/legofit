#ifndef MSPREADER_H
#define MSPREADER_H

#include "typedefs.h"
#include <stdio.h>

// constructor
MsprimeReader *MsprimeReader_new(const char *fname);

// destructor
void MsprimeReader_free(MsprimeReader *self);

// Rewind input and reset chr and nucpos. Doesn't work
// if input is stdin.
int MsprimeReader_rewind(MsprimeReader *self);

// Move MsprimeReader to next nucleotide site.
int MsprimeReader_next(MsprimeReader *self);

// Return current chromosome.
unsigned MsprimeReader_chr(MsprimeReader *self);

// Return the dimension of the array of samples
int MsprimeReader_sampleDim(MsprimeReader *self);

// Return number of samples from population i.
int MsprimeReader_nsamples(MsprimeReader *self, int i);

// Return frequency of derived allele in sample from population i.
double MsprimeReader_daf(MsprimeReader *self, int i);

// Return pointer to label of population i
const char *MsprimeReader_lbl(MsprimeReader *self, int i);
#endif
