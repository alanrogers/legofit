/**
 * @file tabpat.c
 * @brief Tabulate site pattern frequencies from vcf files.
 */

#include "typedefs.h"
#include "misc.h"
#include "binary.h"
#include "dafreader.h"
#include "strint.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct Stack Stack;

struct Stack {
    int dim, nused;
    tipId_t *buff;  // not locally owned
};

static void usage(void);
static Stack *Stack_new(int dim, tipId_t buff[dim]);
static void Stack_free(Stack *stk);
static void Stack_push(Stack *self, tipId_t x);
static void generatePatterns(int bit,  int npops, Stack *stk, tipId_t pat);
static void parseChromosomeLbls(const char *arg, StrInt *strint);

const char *useMsg =
    "\nUsage: tabpat <chromosomes> <x>=<in1> <y>=<in2> ...\n"
    "   where <chromosomes> is an ordered, comma-separated list\n"
    "   of chromosome labels, <x> and <y> are arbitrary labels, and\n"
    "   <in1> and <in2> are input files in daf format. Writes to standard\n"
    "   output. Labels may not include the character \":\".\n";

static void usage(void) {
    fputs(useMsg, stderr);
    fprintf(stderr,"   Maximum number of input files: %lu.\n",
            8*sizeof(tipId_t));
    exit(1);
}

/// This stack is local to this file. It provides a bounds-controlled
/// interface to an external array, which is passed as an argument, buff,
/// to Stack_new.
static Stack *Stack_new(int dim, tipId_t buff[dim]) {
    Stack *self = malloc(sizeof(Stack));
    CHECKMEM(self);
    self->dim = dim;
    self->buff = buff;
    self->nused = 0;
    return self;
}

/// Frees the stack but not the underlying buffer.
static void Stack_free(Stack *stk) {
    free(stk);
}

/// Add an entry to the stack, checking bounds.
static void Stack_push(Stack *self, tipId_t x) {
    if(self->nused == self->dim) {
        fprintf(stderr,"%s:%s:%d buffer overflow\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }
    self->buff[self->nused++] = x;
}

/// Call as generatePatterns(0, npops, stk, 0);
/// Recursive function, which generates all legal site patterns
/// and pushes them onto a stack.
static void generatePatterns(int bit, int npops, Stack *stk, tipId_t pat) {
    assert(sizeof(tipId_t) < sizeof (unsigned long long));
    if(bit == npops) {
        // Exclude patterns with 1 bit on, all bits on, or all bits off.
        if(pat!=0                          // exclude if all bits off
		   && !isPow2(pat)                 // exclude if only 1 bit on
		   && pat != (1ULL << npops) -1ULL // exclude if all bits on
			)
            Stack_push(stk, pat);
        return;
    }
    tipId_t on = 1UL << bit;
    generatePatterns(bit+1, npops, stk, pat|on); // curr bit on
    generatePatterns(bit+1, npops, stk, pat);    // curr bit off
}

/// On input, "arg" is a comma-separated list of chromosome labels.
/// On return, the first label is associated (in strint) with index 0,
/// the second with index 1, and so on. If one of the labels is of form
/// x-y, where x and y are integers, then x-y is expanded into a
/// sequence of integers ranging from x to y inclusive. These integers
/// are converted into strings, and then associated with consequtive
/// integers in the usual way.
///
/// Example: if "arg" is "x,2-4,mtdna", then we end up with the
/// following mapping:
///     x -> 0
///     2 -> 1
///     3 -> 2
///     4 -> 3
/// mtdna -> 4
static void parseChromosomeLbls(const char *arg, StrInt *strint) {
    int i, j;
    char chrs[100], *ptr, *token;
    snprintf(chrs, sizeof chrs, "%s", arg);
    if(NULL != strchr(chrs, '=')) {
        fprintf(stderr,"Bad list of chromosomes: %s\n", chrs);
        usage();
    }
    i = 0;
    ptr = chrs;
    fprintf(stderr,"Expecting chromosomes:");
    while( (token = strsep(&ptr, ",")) != NULL) {
        char *dash = strchr(token, '-');
        int isinteger=1;
        if(dash) {
            if(dash == token) {
                isinteger = 0;
            }else {
                char *p;
                for(p = token; p < dash; ++p) {
                    if(!isdigit(*p)) {
                        isinteger = 0;
                        break;
                    }
                }
                for(p = dash+1; isinteger && *p != '\0'; ++p) {
                    if(!isdigit(*p)) {
                        isinteger = 0;
                        break;
                    }
                }
            }
        }
        if(dash && isinteger) {  // parse range
            *dash = '\0';
            int from = strtol(token, NULL, 10);
            int to = strtol(dash+1, NULL, 10);
            for(j=from; j <= to; ++j) {
                char curr[10];
                snprintf(curr, sizeof curr, "%d", j);
                fprintf(stderr, " %s", curr);
                StrInt_insert(strint, curr, i++);
            }
        }else{                  // token isn't a range
            fprintf(stderr, " %s", token);
            StrInt_insert(strint, token, i++);
        }
    }
    putc('\n', stderr);
}

int main(int argc, char **argv) {
    int i, j;
    int n = argc-2; // number of input files
    char *poplbl[n];
    char *fname[n];
    LblNdx lndx;
    LblNdx_init(&lndx);
	DAFReader *r[n];

    if(n == 0)
        usage();

    // Number of inputs can't exceed number of bits in an object of
    // type tipId_t.
    if(n > 8*sizeof(tipId_t)) {
        fprintf(stderr,"Error: %d input files. Max is %lu.\n",
                n, 8*sizeof(tipId_t));
        usage();
    }

    // 1st argument is comma-separated list of chromosome labels
    StrInt *strint = StrInt_new();
    parseChromosomeLbls(argv[1], strint);

    // Parse remaining arguments, each of which should be of form
    // x=foo, where x is an arbitrary label and foo is the name of an
    // input file.
    for(i=0; i<n; ++i) {
        fname[i] = poplbl[i] = argv[i+2];
        (void) strsep(fname+i, "=");
        if(fname[i] == NULL
           || poplbl[i] == NULL
           || strlen(poplbl[i])==0
           || strlen(fname[i])==0
           || strchr(poplbl[i], ':') != NULL)
            usage();
        LblNdx_addSamples(&lndx, 1, poplbl[i]);
		r[i] = DAFReader_new(fname[i]);
    }

    printf("# Population labels:\n");
    for(i=0; i<n; ++i)
        printf("# %4s = %s\n", poplbl[i], fname[i]);

	// make sure labels are all different
	for(i=1; i<n; ++i)
		for(j=0; j<i; ++j)
			if(0 == strcmp(poplbl[i], poplbl[j])) {
				fprintf(stderr,"Error: duplicate labels on command line.\n");
				fprintf(stderr,"       duplicated label: %s\n", poplbl[i]);
				exit(EXIT_FAILURE);
			}

    unsigned long npat = (1UL<<n) - n - 2; // number of site patterns
    printf("# Number of site patterns: %lu\n", npat);
    tipId_t pat[npat];
	double  patCount[npat];
    int lblsize = 100;
    char lblbuff[lblsize];
	memset(patCount, 0, sizeof(patCount));

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

	fflush(stdout);
	unsigned long nsnps = 0;
	// Iterate through daf files
	while(EOF != DAFReader_multiNext(n, r, strint)) {

        // Skip loci at which data sets disagree about which allele
        // is derived and which ancestral.
        if(!DAFReader_allelesMatch(n, r))
            continue;

		// p and q are frequencies of derived and ancestral alleles
		double p[n], q[n];
		for(j=0; j < n; ++j) {
			p[j] = DAFReader_daf(r[j]);  // derived allele freq
			q[j] = 1-p[j];
		}

		// Contribution of current snp to each site pattern.  Inner
		// loop considers each bit in current pattern.  If that bit is
		// on, multiply z by the derived allele frequency, p. If
		// that bit is off, multiply by q=1-p. In the end, z is Prod
		// p[j]^bit[j] * q[j]^(1-bit[j]) where bit[j] is the value (0
		// or 1) of the j'th bit.
		for(i=0; i < npat; ++i) {
			tipId_t pattern = pat[i];
			double z = 1.0;
			for(j=0; j < n; ++j) {
				if(pattern & 1u)
					z *= p[j];
				else
					z *= q[j];
				pattern >>= 1u;
			}
			assert( 0 == (pattern&1) );
			patCount[i] += z;
		}
		++nsnps;
	}
	printf("# Tabulated %lu SNPs\n", nsnps);

    // print labels and binary representation of site patterns
	printf("# %13s %20s\n", "SitePat", "E[count]");
    for(i=0; i<npat; ++i) {
        printf("%15s %20.7lf\n",
			   patLbl(lblsize, lblbuff,  pat[i], &lndx),
			   patCount[i]);
    }

    for(i=0; i<n; ++i)
		DAFReader_free(r[i]);
    StrInt_free(strint);
    return 0;
}

