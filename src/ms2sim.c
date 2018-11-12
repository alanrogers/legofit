/**
 * @file ms2sim.c
 * @page ms2sim
 * @brief ms2sim translates output from Dick Hudson's `ms` program
 * into "sim" format.
 *
 * # ms2sim, a program that converts `ms` output into sim format
 *
 * This program is written in C, because the conversion involves
 * transposing a large matrix, and this is faster in C.  
 *
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "misc.h"
#include "tokenizer.h"
#include "linereader.h"
#include "uintqueue.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>

typedef struct MSSample MSSample;

struct MSSample {
    double age;
    int pop;      // population
    int nsamp;    // number of haploid samples
    MSSample *next;
};

void      usage(void);
unsigned *countSamples(Tokenizer *tkz, int *sampleDim);
MSSample *MSSample_new(double age, int pop, int nsamp);
int       MSSample_compare(MSSample *left, MSSample *right);
MSSample *MSSample_insert(MSSample *node, MSSample *new);
MSSample *MSSample_collapse(MSSample *node);
int       MSSample_length(MSSample *self);
void      MSSample_free(MSSample *self);
void      MSSample_print(MSSample *self, FILE *fp);

MSSample *MSSample_new(double age, int pop, int nsamp) {
    MSSample *self = malloc(sizeof(MSSample));
    CHECKMEM(self);
    memset(self, 0, sizeof(MSSample));
    self->age = age;
    self->pop = pop;
    self->nsamp = nsamp;
    return self;
}

void MSSample_free(MSSample *self) {
    if(self == NULL)
        return;
    MSSample_free(self->next);
    free(self);
}

// Sort by age.
int MSSample_compare(MSSample *left, MSSample *right) {
    if(left->age > right->age)
        return 1;
    if(left->age < right->age)
        return -1;
    // ages equal
    return 0;
}

int MSSample_length(MSSample *self) {
    if(self == NULL)
        return 0;
    return 1 + MSSample_length(self->next);
}

/// Sort by increasing age. Within age categories, nodes that are
/// inserted later sort as though they were older.
MSSample *MSSample_insert(MSSample *node, MSSample *new) {
    if(node == NULL)
        return new;
    int cmp = MSSample_compare(new, node);
    if(cmp >= 0) {
        node->next = MSSample_insert(node->next, new);
        return node;
    }
    new->next = node;
    return new;
}

/// Collapse adjacent nodes with the same age and population.
MSSample *MSSample_collapse(MSSample *node) {
    if(node == NULL)
        return NULL;
    while(node->next
          && node->age == node->next->age
          && node->pop == node->next->pop) {
        node->nsamp += node->next->nsamp;
        MSSample *tmp = node->next;
        node->next = node->next->next;
        free(tmp);
    }
    node->next = MSSample_collapse(node->next);
    if(node->nsamp == 0) {
        MSSample *tmp = node;
        node = node->next;
        free(tmp);
    }
    return node;
}

void MSSample_print(MSSample *self, FILE *fp) {
    if(self == NULL) 
        return;
    fprintf(stderr,"MSSample: age=%0.8lg, pop=%d, nsamp=%d\n",
            self->age, self->pop, self->nsamp);
    MSSample_print(self->next, fp);
}

/// Parse output of ms
void usage(void) {
    fprintf(stderr, "usage: ms2sim\n");
    fprintf(stderr, "Reads from stdin; writes to stdout.\n");

    exit(EXIT_FAILURE);
}

// On input, tkz should point to a tokenized string representing the
// ms command line, and nsamples should point to an int.
//
// The function returns a newly-allocated array of ints, whose dimension
// is *sampleDim, the number of populations specified on the ms command line.
// The i'th entry in this array is the haploid sample size of population i.
//
// On error, the function returns NULL.
unsigned *countSamples(Tokenizer *tkz, int *sampleDim) {
    /*
      Here is the summary, from Hudson's msdoc.pdf, of the order in
      which haplotypes are printed in ms output: "The ancient DNA
      haplotypes are output after the haplotypes sampled in the
      present. nsam parameter is the number of chromosomes sampled at
      the present time. If multiple -eA switches are employed the
      ancient DNA haplotypes are output after the present day
      haplotypes in order from most recent to most ancient. If ancient
      DNA is sampled at the same time from two or more different
      subpopulations, the haplotypes are output in order of the -eA
      switches on the command line."

      In other words, samples are sorted first by age (youngest first)
      and then by the order in which they appear on the command
      line. The modern sample sizes are specified as arguments to -I
      and the ancient ones as arguments to -eA.
    */

    if(strcmp("ms", Tokenizer_token(tkz, 0)) != 0) {
        fprintf(stderr,"%s:%d: input file is not ms output\n",
                __FILE__,__LINE__);
        return NULL;
    }

    int i, j, npops=0, h;
    char *token, *end;
    MSSample *mss = NULL;
    int ntokens = Tokenizer_ntokens(tkz);

    // Read through tokens, looking for -I and -eA. Use these arguments
    // to set npops and push sample sizes into mss.
    i = 1;
    while(i < ntokens) {
        token = Tokenizer_token(tkz, i);
        if(strcmp("-I", token) == 0) {
            assert(npops == 0);
            assert(mss == NULL);
            token = Tokenizer_token(tkz, i+1);
            npops = strtol(token, &end, 10);
            if(end==token || npops <= 0) {
                fprintf(stderr,"%s:%d: illegal number of populations: %s\n",
                        __FILE__,__LINE__, token);
                goto fail;
            }
            for(j=0; j<npops; ++j) {
                token = Tokenizer_token(tkz, i+j+2);
                h = strtol(token, &end, 10);
                if(end==token || h < 0) {
                    fprintf(stderr,"%s:%d: illegal number of populations: %s\n",
                            __FILE__,__LINE__, token);
                    goto fail;
                }
                // time=0.0, pop=j+1, sample size = h
                MSSample *new = MSSample_new(0.0, j+1, h);
                mss = MSSample_insert(mss, new);
            }
            i += npops + 2;
            continue;
        }else if(strcmp("-eA", token) == 0) {
            token = Tokenizer_token(tkz, i+1);
            double t = strtod(token, &end);
            if(end==token || t < 0.0) {
                fprintf(stderr,"%s:%d: illegal time: %s\n",
                        __FILE__,__LINE__, token);
                goto fail;
            }
            token = Tokenizer_token(tkz, i+2);
            int pop = strtol(token, &end, 10);
            if(end==token || pop < 1) {
                fprintf(stderr,"%s:%d: illegal population: %s\n",
                        __FILE__,__LINE__, token);
                goto fail;
            }
            token = Tokenizer_token(tkz, i+3);
            h = strtol(token, &end, 10);
            if(end==token || h <= 0) {
                fprintf(stderr,"%s:%d: illegal sample size: %s\n",
                        __FILE__,__LINE__, token);
                goto fail;
            }
            // time=t, pop=pop, sample size = h
            MSSample *new = MSSample_new(t, pop, h);
            mss = MSSample_insert(mss, new);
            i += 4;
            continue;
        }
        i += 1;
    }

    mss = MSSample_collapse(mss);
    fprintf(stderr,"Order of samples in output:\n");
    MSSample_print(mss, stderr);
    npops = MSSample_length(mss);

    *sampleDim = npops;
    unsigned *nsamples = malloc(npops * sizeof(nsamples[0]));
    CHECKMEM(nsamples);
    MSSample *node = mss;

    for(i=0; i < npops; ++i) {
        assert(node);
        nsamples[i] = node->nsamp;
        node = node->next;
    }
    MSSample_free(mss);

    return nsamples;

 fail:
    MSSample_free(mss);
    return NULL;
}

int main(int argc, char **argv) {
    size_t      buffsize = 1<<9;
    long        maxtokens = 1<<8;
    long        i, j, ntokens;
    long        site, start, pop, seq;
    long        nsite, nseq;
    time_t      currtime = time(NULL);

    printf("# ms2sim was run %s", ctime(&currtime));
    printf("# cmd lines");
    for(i = 0; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');

    /* command-line arguments */
    int n = 0;
    for(i=1; i < argc; ++i) {
        if(argv[i][0] == '-')
            usage();
        ++n;
    }
    if(n == 0)
        usage();

    char       *poplbl[n];  // array of population labels
    LblNdx      lndx;
    LblNdx_init(&lndx);

    // Number of inputs can't exceed number of bits in an object of
    // type tipId_t.
    if(n > 8 * sizeof(tipId_t)) {
        fprintf(stderr, "Error: %d populations. Max is %lu.\n",
                n, 8 * sizeof(tipId_t));
        usage();
    }

    // Set poplbl array
    for(i=1, j=0; i < argc; ++i) {
        poplbl[j] = argv[i];
        if(strchr(poplbl[j], ':') != NULL) {
            fprintf(stderr,"Label (%s) is illegal because it includes ':'\n",
                    poplbl[j]);
            exit(EXIT_FAILURE);
        }
        LblNdx_addSamples(&lndx, 1, poplbl[j]);
        ++j;
    }

    // allocate memory
    LineReader *lr = LineReader_new(buffsize);
    CHECKMEM(lr);
    Tokenizer *tkz = Tokenizer_new(maxtokens);
    CHECKMEM(tkz);

    char *line;  // not locally owned. Freed by LineReader_free.
    line = LineReader_next(lr, stdin);
    if(line == NULL)
        DIE("No input");
    Tokenizer_split(tkz, line, " \n");
    ntokens = Tokenizer_strip(tkz, " \t\n");

    // First line of input has MS command line
    int sampleDim=0;
    unsigned *nsamples = countSamples(tkz, &sampleDim);

    if(n != sampleDim) {
        fprintf(stderr, 
                "%s:%d: number of labels on command line doesn't match\n"
                " number of subdivisions in ms input.\n",
                __FILE__,__LINE__);
        fprintf(stderr,"  %d labels != %d subdivisions\n",
                n, sampleDim);
        exit(EXIT_FAILURE);
    }

    printf("npops = %d\n", n);
    printf("pop sampsize\n");
    for(i=0; i < n; ++i)
        printf("%s %u\n", poplbl[i], nsamples[i]);
    fflush(stdout);

    // Set nseq by counting the samples in each subdivision. This will
    // not equal the 2nd argument of ms if there are any ancient
    // samples.
    nseq = 0;
    for(i = 0; i < n; ++i)
        nseq += nsamples[i];

    // A matrix of dimension nseq X nsite. To access element for
    // given seq and site: m[nsite*seq + site];
    size_t mdim = 0;
    char *m=NULL;

    // Each iteration processes one replicate
    while(1) {
        // set nsite, number of segregating sites in this replicate
        nsite = -1;
        while(nsite == -1) {
            line = LineReader_next(lr, stdin);
            if(line == NULL)
                goto no_more_data;

            Tokenizer_split(tkz, line, " \t");
            ntokens = Tokenizer_strip(tkz, "\n");
            if(ntokens == 0)
                continue;
            if(0 == strcmp("segsites:", Tokenizer_token(tkz, 0))) {
                // got segsites
                nsite = strtol(Tokenizer_token(tkz, 1), NULL, 10);
            }
        }
        if(nsite <= 0) {
            fprintf(stderr, "%s:%d: Illegal number of sites: %ld",
                    __FILE__, __LINE__, nsite);
            exit(EXIT_FAILURE);
        }

        // skip positions line
        line = LineReader_next(lr, stdin);
        if(line == NULL || 0 != strncmp(line, "positions:", 10))
            DIE("positions not found in input");

        size_t need = nseq * nsite * sizeof(m[0]);
        if(need > mdim) {
            mdim = need;
            m = realloc(m, mdim);
            CHECKMEM(m);
        }

        // Read data, putting values into array "m"
        seq = 0;
        while(1) {
            line = LineReader_next(lr, stdin);
            if(line == NULL) {
                if(seq == nseq)
                    goto done_reading_replicate;
                DIE("missing genotypes in input");
            }

            // len is length line excluding trailing whitespace
            long len;
            for(len=0; line[len] != '\0'; ++len)
                ;
            while(isspace(line[len - 1])) {
                --len;
                line[len] = '\0';
            }
            assert(len == strlen(line));

            if(len == 0)
                goto done_reading_replicate;

            if(seq >= nseq) {
                fprintf(stderr,"%s:%d: too many sequences: %ld (max=%ld)\n",
                        __FILE__, __LINE__, seq+1, nseq);
                fprintf(stderr,"line:");
                for(j=0; j<50 && line[j]!='\0'; ++j)
                    putc(line[j], stderr);
                putc('\n', stderr);
                exit(EXIT_FAILURE);
            }

            if(len != nsite) {
                fprintf(stderr,
                        "%s:%d input line has %ld chars; should have %ld\n",
                        __FILE__, __LINE__, len, nsite);
                exit(EXIT_FAILURE);
            }

            memcpy(m + nsite * seq, line, nsite * sizeof(m[0]));
            seq += 1;
        }
    done_reading_replicate:
        assert(seq == nseq);

        // output
        double nderived, daf;
        for(site = 0; site < nsite; ++site) {
            for(pop=start=0; pop < n; ++pop) {
                nderived = 0.0;
                for(seq=start;  seq < start + nsamples[pop]; ++seq) {
                    assert(seq < nseq);
                    if('1' == m[nsite*seq + site])
                        ++nderived;
                }
                daf = nderived / ((double) nsamples[pop]);
                if(pop > 0)
                    putchar(' ');
                printf("%lf", daf);
                start += nsamples[pop];
            }
            putchar('\n');
        }
        fflush(stdout);
    }
 no_more_data:
    if(m != NULL)
        free(m);

    LineReader_free(lr);
    Tokenizer_free(tkz);

    return 0;
}
