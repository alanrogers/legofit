#ifndef RAFREADER_INCLUDED
#define RAFREADER_INCLUDED

#include "typedefs.h"
#include <stdio.h>

#define RAFSTRSIZE 20

struct RAFReader {
    char       *fname;
    FILE       *fp;
    Tokenizer  *tkz;

    // properties of current snp
    long        snpid;          // 0-based index of current snp
    char ref[RAFSTRSIZE];       // reference allele
    char alt[RAFSTRSIZE];       // alternate alleles
    char chr[RAFSTRSIZE];       // chromosome
    unsigned long nucpos;       // nucleotide position from raf file
    double      p;              // frequency of ancestral allele
};

RAFReader  *RAFReader_new(const char *fname);
void        RAFReader_clearChromosomes(int n, RAFReader *r[n]);
void        RAFReader_free(RAFReader * self);
int         RAFReader_next(RAFReader * self);
double      RAFReader_raf(RAFReader *r);
int         RAFReader_allelesMatch(int n, RAFReader *r[n]);
void        RAFReader_printHdr(FILE *fp);
void        RAFReader_print(RAFReader *r, FILE *fp);
int         RAFReader_rewind(RAFReader *self);
int         RAFReader_multiNext(int n, RAFReader *r[n]);
static inline const char *RAFReader_chr(RAFReader *self);
static inline unsigned long RAFReader_nucpos(RAFReader *self);


/// Return const pointer to label of current chromosome.
static inline const char *RAFReader_chr(RAFReader *self) {
	return self->chr;
}

/// Return position of current nucleotide site
static inline unsigned long RAFReader_nucpos(RAFReader *self) {
    return self->nucpos;
}
#endif
