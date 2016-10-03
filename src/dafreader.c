#include "dafreader.h"
#include "tokenizer.h"
#include "misc.h"
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <ctype.h>

#define MAXFIELDS 5

int iscomment(const char *s);

DAFReader *DAFReader_new(const char *fname) {
    DAFReader *self = malloc(sizeof(*self));
    CHECKMEM(self);
    memset(self, 0, sizeof(DAFReader));
    self->fname = strdup(fname);
    CHECKMEM(self->fname);
    self->fp = fopen(self->fname, "r");
    if(self->fp == NULL) {
        fprintf(stderr,"%s:%s:%d: can't open \"%s\" for input.\n",
                __FILE__,__func__,__LINE__, self->fname);
        exit(EXIT_FAILURE);
    }
    self->tkz = Tokenizer_new(MAXFIELDS);
    self->snpid = -1;
	self->p = strtod("NaN", NULL);
    return self;
}

void DAFReader_free(DAFReader *self) {
    fclose(self->fp);
    free(self->fname);
    Tokenizer_free(self->tkz);
    free(self);
}

// Return 1 if first non-white character in string is '#';
// return 0 otherwise.
int iscomment(const char *s) {
    int rval;
    while(*s != '\0' && isspace(*s))
        ++s;
    rval = (*s == '#');
    return rval;
}

// Read the next snp. Return 0 on success; EOF on end of file.
// Aborts with message if other errors occur.
int DAFReader_next(DAFReader *self) {
    int ntokens1;
    int ntokens;
    int status;
    char buff[100];

    // Find a line of input
    while(1){
        if(fgets(buff, sizeof(buff), self->fp) == NULL)
            return EOF;
        if(NULL == strchr(buff, '\n') && !feof(self->fp)) {
            fprintf(stderr, "%s:%d: Buffer overflow. size=%zu\n",
                    __FILE__,__LINE__, sizeof(buff));
            exit(EXIT_FAILURE);
        }
        if(iscomment(buff))
            continue;
        ntokens1 = Tokenizer_split(self->tkz, buff, " ");
        ntokens = Tokenizer_strip(self->tkz, " \n");
        if(ntokens > 0)
            break;
    }

    if(ntokens != 5) {
        fprintf(stderr,"%s:%d: Each line of .daf file must have 5 tokens,"
                " but current line has %d.\n", __FILE__,__LINE__, ntokens);
        fprintf(stderr,"ntokens1=%d\n", ntokens1);
        fprintf(stderr,"buff: %s\n", buff);
        Tokenizer_printSummary(self->tkz, stderr);
        Tokenizer_print(self->tkz, stderr);
        exit(EXIT_FAILURE);
    }

    ++self->snpid;

    // Chromosome
    char prev[DAFSTRSIZE];
    assert(sizeof prev == sizeof self->chr);
    memcpy(prev, self->chr, sizeof prev);
    status = snprintf(self->chr, sizeof self->chr, "%s",
                      Tokenizer_token(self->tkz, 0));
    if(status >= sizeof self->chr) {
        fprintf(stderr,"%s:%d: chromosome name too long: %s\n",
                __FILE__, __LINE__, Tokenizer_token(self->tkz, 0));
        exit(EXIT_FAILURE);
    }
    int diff = strcmp(prev, self->chr);
    if(diff > 0) {
        fprintf(stderr,"%s:%s:%d: Chromosomes missorted in input.\n",
                __FILE__,__func__,__LINE__);
        fprintf(stderr,"          \"%s\" precedes \"%s\".\n",
                prev, self->chr);
        exit(EXIT_FAILURE);
    }

    // Nucleotide position
    self->nucpos = strtoul(Tokenizer_token(self->tkz, 1), NULL, 10);

    // Ancestral allele
    status = snprintf(self->aa, sizeof(self->aa), "%s",
                      Tokenizer_token(self->tkz, 2));
    strlowercase(self->aa);
    if(strlen(self->aa) != 1 || strchr("atgc", *self->aa) == NULL) {
        fprintf(stderr,"%s:%d: Ancestral allele must be a single nucleotide."
                " Got: %s.\n", __FILE__,__LINE__, self->aa);
        exit(EXIT_FAILURE);
    }
            
    // Derived allele
    snprintf(self->da, sizeof(self->da), "%s", Tokenizer_token(self->tkz, 3));
    strlowercase(self->da);
    if(strlen(self->da) != 1 || strchr("atgc", *self->da) == NULL) {
        fprintf(stderr,"%s:%d: Derived allele must be a single nucleotide."
                " Got: %s.\n", __FILE__,__LINE__, self->da);
        exit(EXIT_FAILURE);
    }

    // Derived allele frequency
    self->p = strtod(Tokenizer_token(self->tkz, 4), NULL);

    return 0;
}

// Returns 0 on success, -1 on failure
int DAFReader_rewind(DAFReader *self) {
    return fseek(self->fp, 0L, SEEK_SET);
}

const char *DAFReader_chr(DAFReader *self) {
	return self->chr;
}

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y) ((X) > (Y) ? (Y) : (X))

#define PRCHR do{                                                       \
        int k;                                                          \
        for(k=0; k<n; ++k)                                              \
            fprintf(stderr," %s%s", k==imaxchr ? "*" : "", r[k]->chr);   \
        putc('\n', stderr);                                             \
}while(0)


/// Advance an array of DAFReaders to the next shared position.
/// Return 0 on success or EOF on end of file.
int DAFReader_multiNext(int n, DAFReader *r[n]) {
    int i;
    unsigned long maxnuc=0, minnuc=ULONG_MAX;
	int imaxchr;       // index of reader with maximum chromosome position
    int onSameChr;    // indicates whether all readers are on same chromosome.
    int diff;
    char old[3][20] = {{'\0'},{'\0'},{'\0'}};

	// Find initial min and max position and chromosome.
    if(EOF == DAFReader_next(r[0]))
        return EOF;
    strcpy(old[0], r[0]->chr);
    imaxchr = 0;
    onSameChr=1;
    for(i=1; i<n; ++i) {
        if(EOF == DAFReader_next(r[i]))
            return EOF;
        strcpy(old[i], r[i]->chr);

        diff = strcmp(r[i]->chr, r[imaxchr]->chr);
        if(diff > 0) {
            onSameChr = 0;
            imaxchr=i;
        }else if(diff < 0)
            onSameChr = 0;
    }

	// Loop until both chr and position are homogeneous.
    do {
        // get them all on the same chromosome
        while(!onSameChr) {
            onSameChr=1;
            for(i=0; i<n; ++i) {
                if(i==imaxchr)
                    continue;
                while((diff=strcmp(r[i]->chr, r[imaxchr]->chr)) < 0) {
                    if(EOF == DAFReader_next(r[i]))
                        return EOF;
                    if(0 != strcmp(old[i], r[i]->chr)) {
                        fprintf(stderr,"%d: old[%d]=%s chr=%s imax=%d max=%s\n",
                                __LINE__, i, old[i], r[i]->chr, imaxchr, r[imaxchr]->chr);
                        strcpy(old[i], r[i]->chr);
                        PRCHR;
                    }
                }
                assert(diff >= 0);
                if(diff > 0) {
                    imaxchr = i;
                    onSameChr=0;
                }
            }
        }

        //PRCHR;

        assert(onSameChr);
        maxnuc = minnuc = r[0]->nucpos;
        for(i=1; i<n; ++i) {
            maxnuc = MAX(maxnuc, r[i]->nucpos);
            minnuc = MIN(minnuc, r[i]->nucpos);
        }

		// Now get them all on the same position. Have
		// to keep checking chr in case one file moves
		// to another chromosome.
        for(i=0; onSameChr && i<n; ++i) {
			// Increment each reader so long as we're all on the same
			// chromosome and the reader's nucpos is low.
            while(onSameChr && r[i]->nucpos < maxnuc) {
                if(EOF == DAFReader_next(r[i]))
                    return EOF;
                if(0 != strcmp(old[i], r[i]->chr)) {
                    fprintf(stderr,"%d: old[%d]=%s chr=%s imax=%d max=%s\n",
                            __LINE__, i, old[i], r[i]->chr, imaxchr, r[imaxchr]->chr);
                    strcpy(old[i], r[i]->chr);
                    PRCHR;
                }
                diff = strcmp(r[i]->chr, r[imaxchr]->chr);
                if(diff > 0) {
                    onSameChr = 0;
                    fprintf(stderr,"%d: imaxchr=%d\n",__LINE__,imaxchr);
                    imaxchr=i;
                }else if(diff < 0)
                    onSameChr = 0;
            }
            if(!onSameChr) {
                fprintf(stderr,"Switched chr[%d]: %s; imaxchr=%d\n", i, r[i]->chr, imaxchr);
                PRCHR;
            }
        }

        //PRCHR;
    }while(!onSameChr || minnuc!=maxnuc);

    return 0;
}

// Return 1 if ancestral and derived alleles of all readers match; 0 otherwise
int DAFReader_allelesMatch(int n, DAFReader *r[n]) {
    int aa = *r[0]->aa;
    int da = *r[0]->da;
    int i;
    for(i=1; i<n; ++i) {
        if(aa != *r[i]->aa || da != *r[i]->da)
            return 0;
    }
    return 1;
}

void DAFReader_printHdr(FILE *fp) {
    fprintf(fp, "%50s %5s %10s %2s %2s %8s\n",
            "file", "chr", "pos", "aa", "da", "daf");
}

void DAFReader_print(DAFReader *r, FILE *fp) {
    assert(r->fname);
    fprintf(fp,"%50s %5s %10lu %2s %2s %8.6lf\n",
            r->fname, r->chr, r->nucpos, r->aa, r->da, r->p);
}

double      DAFReader_daf(DAFReader *r) {
	assert(r->p >= 0.0 && r->p <= 1.0);
	return r->p;
}
