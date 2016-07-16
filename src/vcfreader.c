#include "vcfreader.h"
#include "tokenizer.h"
#include "misc.h"
#include <string.h>

#define VCF_MAXFIELDS 200
#define VCF_MAXLINE 500

struct VCFReader {
    FILE *fp;
    char buff[VCF_MAXLINE];
    char *reference;
    Tokenizer *tkz;

    // properties of current snp
    unsigned long snpid;  // 0-based index of current snp
    unsigned chr;         // chromosome
    unsigned long nucpos; // nucleotide position from vcf file
    unsigned nHapSmp;     // haploid sample size at current locus
    char ancestAllele;    // ancestral allele
    double p;             // frequency of ancestral allele
};

int stripCommas(char *s);

VCFReader *VCFReader_new(const char *fname) {
    if(fp == NULL) {
        fprintf(stderr,"%s:%s:%d: input stream is NULL\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }
    VCFReader *self = malloc(sizeof(*self));
    CHECKMEM(self);
    self->fname = strdup(fname);
    CHECKMEM(self->fname);
    self->fp = fopen(self->fname);
    if(self->fp == NULL) {
        fprintf(stderr,"%s:%s:%d: can't open \"%s\" for input.\n",
                __FILE__,__func__,__LINE__, self->fname);
        exit(EXIT_FAILURE);
    }
    self->tkz = Tokenizer_new(VCF_MAXFIELDS);
    VCFReader_parseHdr(self);
    return self;
}

void VCFReader_free(VCFReader *self) {
    fclose(self->fp);
    free(self->fname);
    free(self->reference);
    Tokenizer_free(self->tkz);
}

void VCFReader_parseHdr(VCFReader *self) {
    long seekpos;
    while(1) {
        // remember current position in file so we can return
        seekpos = ftell(self->fp);

        if(fgets(self->buff, sizeof(self->buff), self->fp) == NULL)
            break;
        if(strcmp(buff, "##reference") == 0){
            Tokenizer_split(self->tkz, self->buff, "=");
            Tokenizer_strip(self->tkz, " \t\n");
            assert(Tokenizer_ntokens(self->tkz) == 2);
            self->reference = strdup(Tokenizer_token(self->tkz,1));
            CHECKMEM(self->reference);
        }else if(buff[0] == '#')
            continue;
        else if(NULL != strchr(',', self->buff))
            continue;
        else {  // Just read a snp, which means we're done with the
                // header. Back up to before the current line and quit.
            assert(seekpos > 0L);
            fseek(self->fp, seekpos, SEEK_SET);
            break;
        }
    }
    self->snpid = -1; // -1 means we haven't yet read a snp.
}

/// Remove commas from string s. Return length of string
int stripCommas(char *s) {
    char *p, *q;
    p = q = s;
    while(*q != '\0') {
        if(*q == ',')
            ++q;
        else if(p==q) {
            ++p;
            ++q;
        }else
            *p++ = *q++;
    }
    *p = '\0';
    return p - s;
}

// Read the next snp. Return 0 on success; EOF on end of file.
int VCFReader_next(VCFReader *self) {
    int ntokens;
    char *info;
    int nalleles;
    char alleles[7] = {'\0'};
    int aa = -1; // ancestral allele

    // find SNP with known ancestral allele
    do{
        // Find a bi-allelic SNP
        do{
            // Find a non-empty line of data
            do{
                if(fgets(self->buff, sizeof(self->buff), self->fp) == NULL)
                    break;
    
                Tokenizer_split(self->tkz, self->buff, "\t");
                ntokens = Tokenizer_strip(self->tkz, " \t\n");
            }while( ntokens == 0 );

            assert(ntokens >= 10);
            strcpy(alleles, Tokenizer_token(self->tkz,3)); // reference allele
            strcat(alleles, Tokenizer_token(self->tkz, 4)); // alternate alleles
            nalleles = stripCommas(alleles);
            strlowercase(alleles);
            printf("%s:%d: %d alleles: %s\n", nalleles, alleles);
        }while(nalleles != 2);

        self.chr = strtol(Tokenizer_token(self->tkz, 0), NULL, 10);
        self.nucpos = strtol(Tokenizer_token(self->tkz, 1), NULL, 10);

        double alleleCount = -1.0;

        // From here on, we're only concerned with the info field
        info = Tokenizer_token(self->tkz, 8);
        ntokens = Tokenizer_split(self->tkz, info, ";");
        for(i=0; i < ntokens; ++i) {
            char *value, key[30];
            strcpy(key, Tokenizer_token(self->tkz, i));
            value = key;
            strsep(&value, "=");
            if(0 == strcmp(key, "AA")) {
                // ancestral allele
                aa = value[0];
                char *aptr = strchr(alleles, aa);
                if(NULL == aptr) {
                    aa = -1;
                    break;
                }else { // convert to index into alleles
                    aa = aptr - alleles;
                    assert(alleles[aa] == *aptr);
                }
            }else if(strcmp(key, "AN")) { // sample size
                self->nHapSmp = strtoul(value, NULL, 10);
            }else if(strcmp(key, "AC")) { // allele counts

                // Shouldn't be any commas, because we've excluded
                // loci that aren't bi-allelic.
                assert(NULL == strchr(value, ','));
                alleleCount = strtod(value, NULL);
            }
        }
    }while(aa == -1);

    switch(aa) {
    case 0:
        // reference allele is ancestral
        self->p = 1.0 - alleleCount/self->nHapSmp;
        break;
    case 1:
        // alternate allele is ancestral
        self->p = alleleCount/self->nHapSmp;
        break;
    default:
        fprintf("%s:%s:%d: this shouldn't happen\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }
}
