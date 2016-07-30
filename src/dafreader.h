#ifndef DAFREADER_INCLUDED
#define DAFREADER_INCLUDED

#include "typedefs.h"
#include <stdio.h>

struct DAFReader {
    char       *fname;
    FILE       *fp;
    Tokenizer  *tkz;

    // properties of current snp
    long        snpid;          // 0-based index of current snp
    char aa[10], da[10];        // ancestral and derived alleles
    unsigned    chr;            // chromosome
    unsigned long nucpos;       // nucleotide position from daf file
    double      p;              // frequency of ancestral allele
};

DAFReader  *DAFReader_new(const char *fname);
void        DAFReader_free(DAFReader * self);
int         DAFReader_next(DAFReader * self);
int         DAFReader_multiNext(int n, DAFReader * r[n]);
double      DAFReader_daf(DAFReader *r);
int         DAFReader_allelesMatch(int n, DAFReader *r[n]);
void        DAFReader_printHdr(FILE *fp);
void        DAFReader_print(DAFReader *r, FILE *fp);

#endif
