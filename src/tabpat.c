/**
 * @file tabpat.c
 * @brief Tabulate site pattern frequencies from .daf files.
 * @copyright Copyright (c) 2016 Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "binary.h"
#include "boot.h"
#include "dafreader.h"
#include "misc.h"
#include "strint.h"
#include "typedefs.h"
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXCHR 24  // maximum number of chromosomes

typedef struct Stack Stack;

struct Stack {
    int dim, nused;
    tipId_t *buff;  // not locally owned
};

static void usage(void);
static Stack *Stack_new(int dim, tipId_t buff[dim]);
static void Stack_free(Stack *stk);
static void Stack_push(Stack *self, tipId_t x);
static void generatePatterns(int bit,  int npops, Stack *stk, tipId_t pat,
                             int doSing);
static void parseChromosomeLbls(const char *arg, StrInt *strint);

const char *useMsg =
    "\nUsage: tabpat [options] <x>=<in1> <y>=<in2> ...\n"
    "   where <x> and <y> are arbitrary labels, and <in1> and <in2> are input\n"
	"   files in daf format. Writes to standard output."
	" Labels may not include\n"
	"   the character \":\".";

static void usage(void) {
    fputs(useMsg, stderr);
    fprintf(stderr," Maximum number of input files: %lu.\n",
            8*sizeof(tipId_t));
	fputs("\nOptions may include:\n", stderr);
	tellopt("-f <name> or --bootfile <name>",
			"Bootstrap output file basename. Def: legoboot.");
	tellopt("-r <x> or --bootreps <x>",
			"# of bootstrap replicates. Def: 0");
	tellopt("-b <x> or --blocksize <x>",
			"# of SNPs per block in moving-blocks bootstrap. Def: 0.");
	tellopt("-1 or --singletons", "Use singleton site patterns");
    tellopt("-m or --logMismatch", "log AA/DA mismatches to tabpat.log");
    tellopt("-F or --logFixed", "log fixed sites to tabpat.log");
    tellopt("-h or --help", "Print this message");
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
        fprintf(stderr,"%s:%s:%d ERR: buffer overflow\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }
    self->buff[self->nused++] = x;
}

/// Call as generatePatterns(0, npops, stk, 0);
/// Recursive function, which generates all legal site patterns
/// and pushes them onto a stack.
static void generatePatterns(int bit, int npops, Stack *stk, tipId_t pat,
                             int doSing) {
    assert(sizeof(tipId_t) < sizeof (unsigned long long));
    if(bit == npops) {
        // Recursion stops here. If current pattern is
        // legal, then push it onto the stack. Then return.
        
        // Exclude patterns with all bits on, or all bits off.
        if(pat==0 || pat == (1ULL << npops) -1ULL)
            return;
        // Exclude singleton patterns unless "doSing" is true.
        if(!doSing && isPow2(pat))
            return;
        Stack_push(stk, pat);
        return;
    }
    tipId_t on = 1UL << bit;
    generatePatterns(bit+1, npops, stk, pat|on, doSing); // curr bit on
    generatePatterns(bit+1, npops, stk, pat, doSing);    // curr bit off
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
			int inc = (to < from ? -1 : 1);
			to += inc; // 1 past final position
            for(j=from; j != to; j += inc) {
                char curr[10];
                snprintf(curr, sizeof curr, "%d", j);
                fprintf(stderr, " %s", curr);
				errno = 0;
                StrInt_insert(strint, curr, i++);
				if(errno) {
					fprintf(stderr,"\n%s:%d: ERR:"
                            " Duplicate chromosome label: %s\n",
							__FILE__,__LINE__, curr);
					exit(EXIT_FAILURE);
				}
            }
        }else{                  // token isn't a range
            fprintf(stderr, " %s", token);
			errno = 0;
            StrInt_insert(strint, token, i++);
			if(errno) {
				fprintf(stderr,"\n%s:%d: ERR: Duplicate chromosome label: %s\n",
						__FILE__,__LINE__, token);
				exit(EXIT_FAILURE);
			}
        }
    }
    putc('\n', stderr);
}

int main(int argc, char **argv) {
    int i, j, status, optndx;
    int doSing=0;       // nonzero means use singleton site patterns
    long bootreps = 0;
    long blocksize = 300;
    StrInt *strint = StrInt_new();
    char bootfname[FILENAMESIZE] = { '\0' };
    const char *logfname = "tabpat.log";
    int logMismatch=0, logFixed=0;
    FILE *logfile = NULL;

    static struct option myopts[] = {
        // {char *name, int has_arg, int *flag, int val}
        {"bootfile", required_argument, 0, 'f'},
        {"bootreps", required_argument, 0, 'r'},
        {"blocksize", required_argument, 0, 'b'},
        {"singletons", no_argument, 0, '1'},
        {"logMismatch", no_argument, 0, 'm'},
        {"logFixed", no_argument, 0, 'F'},
        {"help", no_argument, 0, 'h'},
//        {"threads", required_argument, 0, 't'},
        {NULL, 0, NULL, 0}
    };

    // command line arguments
    for(;;) {
        i = getopt_long(argc, argv, "b:c:f:hr:t:mFv1", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            usage();
            break;
        case 'b':
            blocksize = strtod(optarg, NULL);
			if(blocksize <= 0) {
				fprintf(stderr,
						"%s:%d: bad argument to -b or --blocksize: \"%s\"\n",
						__FILE__,__LINE__,optarg);
				usage();
			}
            break;
		case 'f':
			status = snprintf(bootfname, sizeof bootfname, "%s", optarg);
			if(status >= sizeof bootfname) {
				fprintf(stderr,"%s:%d: ERR: Filename %s is too large."
						" Max: %zu\n",
						__FILE__,__LINE__, optarg,
						sizeof(bootfname)-1);
				exit(EXIT_FAILURE);
			}
			break;
        case 'h':
            usage();
            break;
        case 'r':
            bootreps = strtol(optarg, NULL, 10);
            break;
        case '1':
            doSing=1;
            break;
        case 'm':
            logMismatch=1;
            break;
        case 'F':
            logFixed=1;
            break;
        default:
            usage();
        }
    }

    // remaining options: input files
    int n = argc-optind; // number of input files
    if(n == 0)
        usage();

    char *poplbl[n];
    char *fname[n];
    LblNdx lndx;
    LblNdx_init(&lndx);
	DAFReader *r[n];

    // Number of inputs can't exceed number of bits in an object of
    // type tipId_t.
    if(n > 8*sizeof(tipId_t)) {
        fprintf(stderr,"Error: %d input files. Max is %lu.\n",
                n, 8*sizeof(tipId_t));
        usage();
    }

    // Parse remaining arguments, each of which should be of form
    // x=foo, where x is an arbitrary label and foo is the name of an
    // input file.
    for(i=0; i<n; ++i) {
        fname[i] = poplbl[i] = argv[i+optind];
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

    if(logMismatch || logFixed) {
        logfile = fopen(logfname, "w");
        if(logfile==NULL) {
            fprintf(stderr,"Can't write to file \"%s\".\n",logfname);
            exit(EXIT_FAILURE);
        }
    }

	// Default boot file name
    if(bootfname[0] == '\0') {
		const char *defName = "legoboot";
		status = snprintf(bootfname, sizeof bootfname, "%s", defName);
		if(status >= sizeof bootfname) {
			fprintf(stderr,"%s:%d: ERR: Filename %s is too large."
					" Max: %zu\n",
					__FILE__,__LINE__, defName,
					sizeof(bootfname)-1);
			exit(EXIT_FAILURE);
		}
	}

    printf("# Population labels:\n");
    for(i=0; i<n; ++i)
        printf("# %4s = %s\n", poplbl[i], fname[i]);

	// make sure labels are all different
	for(i=1; i<n; ++i)
		for(j=0; j<i; ++j)
			if(0 == strcmp(poplbl[i], poplbl[j])) {
				fprintf(stderr,"ERR: duplicate labels on command line.\n");
				fprintf(stderr,"     duplicated label: %s\n", poplbl[i]);
				exit(EXIT_FAILURE);
			}

    unsigned long npat = (1UL<<n) - 2; // number of site patterns
    if(!doSing)
        npat -= n;
    printf("# %s singleton site patterns.\n",
           (doSing ? "Including" : "Excluding"));
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
        generatePatterns(0, n, stk, 0, doSing);

        Stack_free(stk);
    }

    // Sort site patterns. Major sort is by number of "on" bits,
    // so that doubleton patterns come first, then tripletons, ets.
    // Secondary sort is by order in which labels are listed
    // on the command line.
    qsort(pat, (size_t) npat, sizeof(pat[0]), compare_tipId);
	fflush(stdout);

    // Used by bootstrap
    Boot *boot = NULL;
	int nchr=0;
    char prev[DAFSTRSIZE], chr[DAFSTRSIZE] = {'\0'};
    long nsnp[MAXCHR];
	memset(nsnp, 0, sizeof nsnp);

    // Read the data to get dimensions: number of chromosomes and
    // number of snps per chromosome. Then use these dimensions to
    // allocate a bootstrap object. 
    if(bootreps > 0) {
        fprintf(stderr,"Doing 1st pass through data to get dimensions...\n");

        // First pass through data sets values of
        // nchr
        // nsnp[i] {i=0..nchr-1}
        while(EOF != DAFReader_multiNext(n, r)) {

            // Skip loci at which data sets disagree about which allele
            // is derived and which ancestral.
            if(!DAFReader_allelesMatch(n, r))
                continue;

            assert(strlen(DAFReader_chr(r[0])) < sizeof prev);
            strcpy(prev, chr);
            strcpy(chr, DAFReader_chr(r[0]));
            int diff = strcmp(prev, chr);
            if(diff != 0) {
                StrInt_insert(strint, chr, nchr);
                nsnp[nchr] = 1;
                ++nchr;
            }else
                ++nsnp[nchr-1];
        }

        for(i=0; i<n; ++i) {
            status = DAFReader_rewind(r[i]);
			if(status) {
				fprintf(stderr,"%s:%d: ERR: can't rewind input stream.\n",
						__FILE__,__LINE__);
				fprintf(stderr,"  If --bootreps > 0, inputs must be"
						" files, not pipes.\n");
				exit(EXIT_FAILURE);
			}
		}

        // Allocate Boot structure
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set(rng, (unsigned long) time(NULL));
        boot = Boot_new(nchr, nsnp, bootreps, npat, blocksize, rng);
        gsl_rng_free(rng);
        CHECKMEM(boot);
    }

	unsigned long nsnps = 0;
    long snpndx = -1;

	// Iterate through daf files
    int chrndx=-1, currChr = INT_MAX;
	while(EOF != DAFReader_multiNext(n, r)) {

        // Skip loci at which data sets disagree about which allele
        // is derived and which ancestral.
        if(!DAFReader_allelesMatch(n, r)) {
            if(logMismatch) {
                assert(logfile);
                fprintf(logfile, "Mismatch\n");
                DAFReader_printHdr(logfile);
                for(i=0; i<n; ++i)
                    DAFReader_print(r[i], logfile);
            }
            continue;
        }

        if(bootreps > 0) {
            // chrndx is index of current chromosome
            errno = 0;
            chrndx = StrInt_get(strint, DAFReader_chr(r[0]));
            if(errno) {
                fprintf(stderr,"%s:%d: ERR: missing index for chromosome: %s\n",
                        __FILE__,__LINE__, DAFReader_chr(r[0]));
                exit(EXIT_FAILURE);
            }
            if(chrndx != currChr) {
                currChr = chrndx;
                snpndx = 0;
            }else
                ++snpndx;

#ifndef NDEBUG
			assert(snpndx < nsnp[chrndx]);
#endif
        }

		// p and q are frequencies of derived and ancestral alleles
		double p[n], q[n], minp, maxp;
        minp=1.0;
        maxp=0.0;
		for(j=0; j < n; ++j) {
			p[j] = DAFReader_daf(r[j]);  // derived allele freq
			q[j] = 1-p[j];
            minp = fmin(minp, p[j]);
            maxp = fmax(maxp, p[j]);
		}
        if(maxp==0.0 || minp==1.0) {
            if(logFixed) {
                assert(logfile);
                fprintf(logfile, "All p=%lf\n", maxp==0.0?0.0:1.0);
                DAFReader_printHdr(logfile);
                for(i=0; i<n; ++i)
                    DAFReader_print(r[i], logfile);
            }
            continue;
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
			if(!isfinite(z)) {
				fprintf(stderr,"%s:%d nonfinite z=%lf\n",
						__FILE__,__LINE__,z);
				fprintf(stderr,"   pattern=%d\n", pat[i]);
				for(j=0; j < n; ++j)
					fprintf(stderr,"   %d: p=%lf q=%lf\n",
							j, p[j], q[j]);
			}
			assert( 0 == (pattern&1) );
			patCount[i] += z;
            if(bootreps > 0) {
                assert(snpndx >= 0);
                assert(chrndx >= 0);
                Boot_add(boot, chrndx, snpndx, i, z);
            }
		}
#ifndef NDEBUG
        if(bootreps > 0)
            Boot_sanityCheck(boot,__FILE__,__LINE__);
#endif
		++nsnps;
		errno = 0;
	}
	printf("# Tabulated %lu SNPs\n", nsnps);

	// boottab[i][j] is the count of the j'th site pattern
	// in the i'th bootstrap replicate.
	double bootvals[bootreps];
	double boottab[bootreps][npat];
	memset(boottab, 0, sizeof boottab);

	if(bootreps > 0) {
		printf("# %s = %s\n", "bootstrap output file", bootfname);
#ifndef NDEBUG
		Boot_sanityCheck(boot,__FILE__,__LINE__);
#endif
        // put site pattern counts into matrix boottab.
		for(i=0; i<bootreps; ++i)
			Boot_aggregate(boot, i, npat, boottab[i]);

        // write an output file for each bootstrap replicate
		for(j=0; j < bootreps; ++j) {
            char buff[FILENAMESIZE+3];
            status = snprintf(buff, sizeof buff, "%s%03d",
                              bootfname, j);
			if(status >= sizeof buff)
                DIE("buffer overflow in snprintf");

            FILE *fp = fopen(buff, "w");
            if(fp == NULL)
                DIE("bad fopen");
            fprintf(fp,"# %13s %20s", "SitePat", "E[count]\n");
            for(i=0; i<npat; ++i) {
                fprintf(fp, "%15s %20.7lf\n",
                       patLbl(lblsize, lblbuff,  pat[i], &lndx),
                       boottab[j][i]);
            }
            fclose(fp);
        }
	}

    // print labels and binary representation of site patterns
	printf("# %13s %20s", "SitePat", "E[count]");
	if(bootreps > 0)
		printf(" %12s %12s", "loBnd", "hiBnd");
	putchar('\n');
    for(i=0; i<npat; ++i) {
        printf("%15s %20.7lf",
			   patLbl(lblsize, lblbuff,  pat[i], &lndx),
			   patCount[i]);
		if(bootreps > 0) {
			double lowBnd, highBnd;
			for(j=0; j < bootreps; ++j)
				bootvals[j] = boottab[j][i];
			confidenceBounds(&lowBnd, &highBnd, 0.95, bootreps, bootvals);
			printf(" %12.7lf %12.7lf", lowBnd, highBnd);
		}
        putchar('\n');
    }

    for(i=0; i<n; ++i)
		DAFReader_free(r[i]);
	if(bootreps > 0)
		Boot_free(boot);
    StrInt_free(strint);
    if(logfile)
        fclose(logfile);
    return 0;
}

