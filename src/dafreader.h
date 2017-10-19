#ifndef DAFREADER_INCLUDED
#define DAFREADER_INCLUDED

#include "typedefs.h"
#include <stdio.h>

#define DAFSTRSIZE 20

struct DAFReader {
    char       *fname;
    FILE       *fp;
    Tokenizer  *tkz;

    // properties of current snp
    long        snpid;          // 0-based index of current snp
    char aa[DAFSTRSIZE];        // ancestral alleles
    char da[DAFSTRSIZE];        // derived alleles
    char chr[DAFSTRSIZE];       // chromosome
    unsigned long nucpos;       // nucleotide position from daf file
    double      p;              // frequency of ancestral allele
};

DAFReader  *DAFReader_new(const char *fname);
void        DAFReader_clearChromosomes(int n, DAFReader *r[n]);
void        DAFReader_free(DAFReader * self);
int         DAFReader_next(DAFReader * self);
double      DAFReader_daf(DAFReader *r);
int         DAFReader_allelesMatch(int n, DAFReader *r[n]);
void        DAFReader_printHdr(FILE *fp);
void        DAFReader_print(DAFReader *r, FILE *fp);
int         DAFReader_rewind(DAFReader *self);
int         DAFReader_multiNext(int n, DAFReader *r[n]);
static inline const char *DAFReader_chr(DAFReader *self);
static inline unsigned long DAFReader_nucpos(DAFReader *self);


/// Return const pointer to label of current chromosome.
static inline const char *DAFReader_chr(DAFReader *self) {
	return self->chr;
}

/// Return position of current nucleotide site
static inline unsigned long DAFReader_nucpos(DAFReader *self) {
    return self->nucpos;
}
#endif
