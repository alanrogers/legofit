#include "vcfreader.h"
#include "tokenizer.h"
#include "misc.h"
#include <string.h>
#include <limits.h>
#include <ctype.h>

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
	self->p = strtod("NaN", NULL);
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
		if(NULL == strchr(self->buff, '\n')) {
			fprintf(stderr, "%s:%d: Buffer overflow. size=%zu\n",
					__FILE__,__LINE__, sizeof(self->buff));
			exit(EXIT_FAILURE);
		}
        if(0 == strncasecmp(self->buff, "##reference", 11)){
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
				if(NULL == strchr(self->buff, '\n')) {
					fprintf(stderr, "%s:%d: Buffer overflow. size=%zu\n",
							__FILE__,__LINE__, sizeof(self->buff));
					exit(EXIT_FAILURE);
				}
    
                Tokenizer_split(self->tkz, self->buff, "\t");
                ntokens = Tokenizer_strip(self->tkz, " \t\n");
            }while( ntokens == 0 );

			if(ntokens < 10) {
				fprintf(stderr,"%s:%s:%d: ERR ntokens=%d < 10\n",
						__FILE__,__func__,__LINE__, ntokens);
				Tokenizer_print(self->tkz, stderr);
				exit(1);
			}
            assert(ntokens >= 10);
            strcpy(alleles, Tokenizer_token(self->tkz,3)); // reference allele
            strcat(alleles, Tokenizer_token(self->tkz, 4)); // alternate alleles
            nalleles = stripCommas(alleles);
            strlowercase(alleles);
        }while(nalleles > 2);

		strcpy(self->alleles, alleles);
        self->chr = strtol(Tokenizer_token(self->tkz, 0), NULL, 10);
        self->nucpos = strtol(Tokenizer_token(self->tkz, 1), NULL, 10);

        alleleCount = -1.0;

        // From here on, we're only concerned with the info field
        info = Tokenizer_token(self->tkz, 7);
        ntokens = Tokenizer_split(self->tkz, info, ";");
        for(i=0; i < ntokens; ++i) {
            char *value, key[30];
            strcpy(key, Tokenizer_token(self->tkz, i));
            value = key;
            strsep(&value, "=");
            if(0 == strcmp(key, "AA")) {
                // ancestral allele
                aa = tolower(value[0]);
                char *aptr = strchr(alleles, aa);
                if(NULL == aptr) {
                    self->ancestAllele = aa = -1;
                    break;
                }else { // convert to index into alleles
                    self->ancestAllele = aa = aptr - alleles;
                    assert(alleles[aa] == *aptr);
                }
            }else if(0 == strcmp(key, "AN")) { // sample size
                self->nHapSmp = strtoul(value, NULL, 10);
            }else if(0 == strcmp(key, "AC")) { // allele counts

                // Shouldn't be any commas, because we've excluded
                // loci that aren't bi-allelic.
                assert(NULL == strchr(value, ','));
                alleleCount = strtod(value, NULL);
            }
        }
    }while(aa == -1);

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

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y) ((X) > (Y) ? (Y) : (X))

/// Advance an array of VCFReaders to the next shared position.
/// Return 0 on success or EOF on end of file.
int VCFReader_multiNext(int n, VCFReader *r[n]) {
    int i;
    unsigned long maxnuc=0, minnuc=ULONG_MAX;
	unsigned maxchr=0, minchr = UINT_MAX;

	// Find initial min and max position and chromosome.
    for(i=0; i<n; ++i) {
        if(EOF == VCFReader_next(r[i]))
            return EOF;
        maxnuc = MAX(maxnuc, r[i]->nucpos);
        minnuc = MIN(minnuc, r[i]->nucpos);
		maxchr = MAX(maxchr, r[i]->chr);
		minchr = MIN(minchr, r[i]->chr);
    }

	// Loop until both chr and position are homogeneous.
    while(minchr!=maxchr || minnuc!=maxnuc) {

		// get them all on the same chromosome
		while(minchr!=maxchr) {
			for(i=0; i<n; ++i) {
				while(r[i]->chr < maxchr) {
					if(EOF == VCFReader_next(r[i]))
						return EOF;
				}
			}
			maxchr=minchr=r[0]->chr;
			for(i=1; i<n; ++i) {
				maxchr = MAX(maxchr, r[i]->chr);
				minchr = MIN(minchr, r[i]->chr);
			}
		}

		// Now get them all on the same position. Have
		// to keep checking chr in case one file moves
		// to another chromosome.
        for(i=0; i<n; ++i) {
			// Increment each reader so long as we're all on the same
			// chromosome and the reader's nucpos is low.
            while(r[i]->chr==maxchr && r[i]->nucpos < maxnuc) {
                if(EOF == VCFReader_next(r[i]))
                    return EOF;
            }
        }

		// Recalculate all max and min values.
        maxnuc=minnuc=r[0]->nucpos;
		maxchr=minchr=r[0]->chr;
        for(i=1; i<n; ++i) {
            maxnuc = MAX(maxnuc, r[i]->nucpos);
            minnuc = MIN(minnuc, r[i]->nucpos);
			maxchr = MAX(maxchr, r[i]->chr);
			minchr = MIN(minchr, r[i]->chr);
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
    fprintf(fp,"  %25s: %s\n", "alleles", r->alleles);
    fprintf(fp,"  %25s: %u\n", "haploid sample size", r->nHapSmp);
    fprintf(fp,"  %25s: %d\n", "ancestral allele", r->ancestAllele);
    fprintf(fp,"  %25s: %lf\n", "ancestral allele freq", r->p);
}

double      VCFReader_aaFreq(VCFReader *r) {
	assert(r->p >= 0.0 && r->p <= 1.0);
	return r->p;
}
