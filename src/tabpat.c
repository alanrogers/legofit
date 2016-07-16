/**
 * @file tabpat.c
 * @brief Tabulate site pattern frequencies from vcf files.
 */

#include "typedefs.h"
#include "misc.h"
#include "binary.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct Stack Stack;

struct Stack {
    int dim, nused;
    tipId_t *buff;  // not locally owned
};

void usage(void);
Stack *Stack_new(int dim, tipId_t buff[dim]);
void Stack_free(Stack *stk);
void Stack_push(Stack *self, tipId_t x);
void generatePatterns(int bit,  int npops, Stack *stk, tipId_t pat);

const char *useMsg =
    "Usage: tabpat <x>=<in1.vcf> <y>=<in2.vcf> ...\n"
    "   where <x> and <y> are arbitrary labels, and <in1> and <in2> are\n"
    "   the names of input files in vcf format. Writes to standard output.\n";

void usage(void) {
    fputs(useMsg, stderr);
    exit(1);
}

Stack *Stack_new(int dim, tipId_t buff[dim]) {
    Stack *self = malloc(sizeof(Stack));
    CHECKMEM(self);
    self->dim = dim;
    self->buff = buff;
    self->nused = 0;
    return self;
}

void Stack_free(Stack *stk) {
    free(stk);
}

void Stack_push(Stack *self, tipId_t x) {
    if(self->nused == self->dim) {
        fprintf(stderr,"%s:%s:%d buffer overflow\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }
    self->buff[self->nused++] = x;
}


/// Call as generatePatterns(0, npops, stk, 0UL);
void generatePatterns(int bit,  int npops, Stack *stk, tipId_t pat) {
    if(bit == npops) {
        // Exclude patterns with 1 bit on, all bits on, or all bits off.
        if(pat==0UL || isPow2(pat))
            return;
        unsigned long long allBitsOn = (1ULL << npops) -1;
        if(pat == allBitsOn)
            return;
        Stack_push(stk, pat);
        return;
    }
    tipId_t on = 1UL << bit;
    generatePatterns(bit+1, npops, stk, pat|on); // curr bit on
    generatePatterns(bit+1, npops, stk, pat);    // curr bit off
}

int main(int argc, char **argv) {
    int i;
    int n = argc-1; // number of inputs
    FILE *in[n];
    char *poplbl[n];
    char *fname[n];
    LblNdx lndx;
    LblNdx_init(&lndx);

    if(n == 0)
        usage();

    // Parse arguments, each of which should be of form
    // x=foo, where x is an arbitrary label and foo is the
    // name of an input file.
    for(i=0; i<n; ++i) {
        fname[i] = poplbl[i] = argv[i+1];
        (void) strsep(fname+i, "=");
        if(fname[i] == NULL || strlen(fname[i])==0)
            usage();
        printf("%4s = %s\n", poplbl[i], fname[i]);
        in[i] = fopen(fname[i], "r");
        if(in[i] == NULL) {
            fprintf(stderr,"Couldn't open \"%s\"\n", fname[i]);
            exit(1);
        }
        LblNdx_addSamples(&lndx, 1, poplbl[i]);
    }

    unsigned long npat = (1UL<<n) - n - 2; // number of site patterns
    printf("npat=%lu\n", npat);fflush(stdout);
    tipId_t pat[npat];

    Stack *stk = Stack_new(npat, pat);
    generatePatterns(0, n, stk, 0);

    qsort(pat, (size_t) npat, sizeof(pat[0]), compare_tipId);

    for(i=0; i<npat; ++i) {
        int lblsize = 100;
        char lblbuff[lblsize];
        printf("%-12s: ", patLbl(lblsize, lblbuff,  pat[i], &lndx));
        printBits(sizeof(tipId_t), pat+i, stdout);
    }

    for(i=0; i<n; ++i)
        fclose(in[i]);
    Stack_free(stk);
    return 0;
}

