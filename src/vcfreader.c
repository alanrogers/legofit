#include "vcfreader.h"
#include "tokenizer.h"
#include "misc.h"
#include <string.h>
#include <limits.h>

#define VCF_MAXFIELDS 200

int stripCommas(char *s);

VCFReader *VCFReader_new(const char *fname) {
    VCFReader *self = malloc(sizeof(*self));
    CHECKMEM(self);
    memset(self, 0, sizeof(VCFReader));
    self->fname = strdup(fname);
    CHECKMEM(self->fname);
    self->fp = fopen(self->fname, "r");
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
    if(self->reference != NULL)
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
        if(strcmp(self->buff, "##reference") == 0){
            Tokenizer_split(self->tkz, self->buff, "=");
            Tokenizer_strip(self->tkz, " \t\n");
            assert(Tokenizer_ntokens(self->tkz) == 2);
            self->reference = strdup(Tokenizer_token(self->tkz,1));
            CHECKMEM(self->reference);
        }else if(self->buff[0] == '#')
            continue;
        else if(NULL != strchr(self->buff, ','))
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
    int i, nalleles;
    char alleles[7] = {'\0'};
    int aa = -1; // ancestral allele
    double alleleCount;

    // find SNP with known ancestral allele
    do{
        // Exclude SNPs with >2 alleles
        do{
            // Find a non-empty line of data
            do{
                if(fgets(self->buff, sizeof(self->buff), self->fp) == NULL)
                    return EOF;
    
                Tokenizer_split(self->tkz, self->buff, "\t");
                ntokens = Tokenizer_strip(self->tkz, " \t\n");
            }while( ntokens == 0 );

            if(ntokens < 10) {
                printf("%s:%d: ntokens=%d\n",__FILE__,__LINE__,ntokens);
                Tokenizer_print(self->tkz, stdout);
            }
            assert(ntokens >= 10);
            strcpy(alleles, Tokenizer_token(self->tkz,3)); // reference allele
            strcat(alleles, Tokenizer_token(self->tkz, 4)); // alternate alleles
            nalleles = stripCommas(alleles);
            strlowercase(alleles);
            printf("%s:%d: %d alleles: %s\n", __FILE__,__LINE__,
                   nalleles, alleles);
        }while(nalleles > 2);

        self->chr = strtol(Tokenizer_token(self->tkz, 0), NULL, 10);
        self->nucpos = strtol(Tokenizer_token(self->tkz, 1), NULL, 10);

        alleleCount = -1.0;

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
                self->ancestAllele = aa = value[0];
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

    assert(alleleCount > 0.0);

    ++self->snpid;

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
        fprintf(stderr, "%s:%s:%d: this shouldn't happen\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }
    return 0;
}



/// Advance an array of VCFReaders to the next shared position.
/// Return 0 on success or EOF on end of file.
int VCFReader_multiNext(int n, VCFReader *r[n]) {
    int i;
    unsigned long max=0, min=ULONG_MAX;

    for(i=0; i<n; ++i) {
        if(EOF == VCFReader_next(r[i]))
            return EOF;
        max = fmax(max, r[i]->nucpos);
        min = fmin(min, r[i]->nucpos);
    }

    while(min != max) {
        for(i=0; i<n; ++i) {
            while(r[i]->nucpos < max) {
                if(EOF == VCFReader_next(r[i]))
                    return EOF;
            }
        }
        max=min=r[0]->nucpos;
        for(i=1; i<n; ++i) {
            max = fmax(max, r[i]->nucpos);
            min = fmin(min, r[i]->nucpos);
        }
    }
    return 0;
}

void VCFReader_print(VCFReader *r, FILE *fp) {
    assert(r->fname);
    fprintf(fp,"VCFReader(%s):\n", r->fname);
    fprintf(fp,"  %25s: %s\n", "reference",
            (r->reference ? r->reference : "NULL"));
    fprintf(fp,"  %25s: %ld\n", "snpid", r->snpid);
    fprintf(fp,"  %25s: %u\n", "chromosome", r->chr);
    fprintf(fp,"  %25s: %lu\n", "nucpos", r->nucpos);
    fprintf(fp,"  %25s: %u\n", "haploid sample size", r->nHapSmp);
    fprintf(fp,"  %25s: %c\n", "ancestral allele", r->ancestAllele);
    fprintf(fp,"  %25s: %lf\n", "ancestral allele freq", r->p);
}
