/**
 * @file tabpat.c
 * @brief Tabulate site pattern frequencies from vcf files.
 */

#include "typedefs.h"
#include "misc.h"
#include "binary.h"
#include "vcfreader.h"
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
    "\nUsage: tabpat <x>=<in1> <y>=<in2> ...\n"
    "   where <x> and <y> are arbitrary labels, and <in1> and <in2> are\n"
    "   the names of input files in vcf format. Writes to standard output.\n"
    "   Labels may not include the character \":\".\n";

void usage(void) {
    fputs(useMsg, stderr);
    fprintf(stderr,"   Maximum number of input files: %lu.\n",
            8*sizeof(tipId_t));
    exit(1);
}

/// This stack is local to this file. It provides a bounds-controlled
/// interface to an external array, which is passed as an argument, buff,
/// to Stack_new.
Stack *Stack_new(int dim, tipId_t buff[dim]) {
    Stack *self = malloc(sizeof(Stack));
    CHECKMEM(self);
    self->dim = dim;
    self->buff = buff;
    self->nused = 0;
    return self;
}

/// Frees the stack but not the underlying buffer.
void Stack_free(Stack *stk) {
    free(stk);
}

/// Add an entry to the stack, checking for buffer overflow.
void Stack_push(Stack *self, tipId_t x) {
    if(self->nused == self->dim) {
        fprintf(stderr,"%s:%s:%d buffer overflow\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }
    self->buff[self->nused++] = x;
}

/// Call as generatePatterns(0, npops, stk, 0UL);
/// Recursive function, which generates all legal site patterns
/// and pushes them onto a stack.
void generatePatterns(int bit,  int npops, Stack *stk, tipId_t pat) {
    assert(sizeof(tipId_t) < sizeof (unsigned long long));
    if(bit == npops) {
        // Exclude patterns with 1 bit on, all bits on, or all bits off.
        if(pat!=0 && !isPow2(pat) && pat != (1ULL << npops) -1ULL)
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
    char *poplbl[n];
    char *fname[n];
    LblNdx lndx;
    LblNdx_init(&lndx);
	VCFReader *r[n];

    if(n == 0)
        usage();

    // Number of inputs can't exceed number of bits in an object of
    // type tipId_t.
    if(n > 8*sizeof(tipId_t)) {
        fprintf(stderr,"Error: %d input files. Max is %lu.\n",
                n, 8*sizeof(tipId_t));
        usage();
    }

    // Parse arguments, each of which should be of form
    // x=foo, where x is an arbitrary label and foo is the
    // name of an input file.
    for(i=0; i<n; ++i) {
        fname[i] = poplbl[i] = argv[i+1];
        (void) strsep(fname+i, "=");
        if(fname[i] == NULL
           || poplbl[i] == NULL
           || strlen(poplbl[i])==0
           || strlen(fname[i])==0
           || strchr(poplbl[i], ':') != NULL)
            usage();
        LblNdx_addSamples(&lndx, 1, poplbl[i]);
		r[i] = VCFReader_new(fname[i]);
    }

    printf("Population labels:\n");
    for(i=0; i<n; ++i)
        printf("%4s = %s\n", poplbl[i], fname[i]);

    unsigned long npat = (1UL<<n) - n - 2; // number of site patterns
    printf("Number of site patterns: %lu\n", npat);fflush(stdout);
    tipId_t pat[npat];

    {
        // Stack is a interface to array "pat".
        Stack *stk = Stack_new(npat, pat);

        // Put site patterns into array "pat".
        generatePatterns(0, n, stk, 0);

        Stack_free(stk);
    }

    // Sort site patterns. Major sort is by number of "on" bits,
    // so that doubleton patterns come first, then tripletons, ets.
    // Secondary sort is by order in which labels are listed
    // on the command line.
    qsort(pat, (size_t) npat, sizeof(pat[0]), compare_tipId);

    // print labels and binary representation of site patterns
    for(i=0; i<npat; ++i) {
        int lblsize = 100;
        char lblbuff[lblsize];
        printf("%15s ", patLbl(lblsize, lblbuff,  pat[i], &lndx));
        printBits(sizeof(tipId_t), pat+i, stdout);
    }

	// Iterate through vcf files
	while(EOF != VCFReader_multiNext(n, r)) {
		fputs("#################\n", stdout);
		for(i=0; i<n; ++i)
			VCFReader_print(r[i], stdout);
	}

    for(i=0; i<n; ++i)
		VCFReader_free(r[i]);
    return 0;
}

