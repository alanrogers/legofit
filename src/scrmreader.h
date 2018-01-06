#ifndef SCRMREADER_H
#define SCRMREADER_H

#include "typedefs.h"
#include <stdio.h>

// constructor
ScrmReader *ScrmReader_new(FILE *fp);

// destructor
void ScrmReader_free(ScrmReader *self);

// Rewind input and reset chr and nucpos. Doesn't work
// if input is stdin.
int ScrmReader_rewind(ScrmReader *self);

// Move ScrmReader to next nucleotide site.
int ScrmReader_next(ScrmReader *self);

// Return current chromosome.
unsigned ScrmReader_chr(ScrmReader *self);

// Return current nucleotide position.
unsigned long ScrmReader_nucpos(ScrmReader *self);

// Return number of populations.
int ScrmReader_npops(ScrmReader *self);

// Return number of samples from population i.
int ScrmReader_nsamples(ScrmReader *self, int i);

// Return frequency of derived allele in sample from population i.
double ScrmReader_daf(ScrmReader *self, int i);
#endif
