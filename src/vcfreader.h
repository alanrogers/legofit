#ifndef VCFREADER_INCLUDED
#define VCFREADER_INCLUDED

#include "typedefs.h"
#include <stdio.h>

#define VCF_MAXLINE 1000

struct VCFReader {
    char       *fname;
    FILE       *fp;
    char        buff[VCF_MAXLINE];
    char       *reference;
    Tokenizer  *tkz;

    // properties of current snp
    long        snpid;          // 0-based index of current snp
	char        alleles[5];
    unsigned    chr;            // chromosome
    unsigned long nucpos;       // nucleotide position from vcf file
    unsigned    nHapSmp;        // haploid sample size at current locus
    int         ancestAllele;   // ancestral allele
    double      p;              // frequency of ancestral allele
};

VCFReader  *VCFReader_new(const char *fname);
void        VCFReader_free(VCFReader * self);
void        VCFReader_parseHdr(VCFReader * self);
int         VCFReader_next(VCFReader * self);
int         VCFReader_multiNext(int n, VCFReader * r[n]);
void        VCFReader_print(VCFReader *r, FILE *fp);
double      VCFReader_aaFreq(VCFReader *r);


#endif
