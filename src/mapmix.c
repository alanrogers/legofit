/**
@file mapmix.c
@page mapmix
@brief Map admixture across the genome

# Mapmix: estimate admixture at each nucleotide site

# mapmix --zero <a.legosim> --one <b.legosim> \
#    <x>=<in_1> <y>=<in_2> ... outgroup=<in_K>

**/

#include "branchtab.h"
#include "error.h"
#include "lblndx.h"
#include "misc.h"
#include "rafreader.h"
#include "strdblqueue.h"
#include "typedefs.h"
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef struct Stack Stack;

/// Treat a vector of tipId_t values as a push-down stack.
struct Stack {
    int         dim, nused;
    tipId_t    *buff;           // not locally owned
};

static void usage(void);
static Stack *Stack_new(int dim, tipId_t buff[dim]);
static void Stack_free(Stack * stk);
static void Stack_push(Stack * self, tipId_t x);
static void generatePatterns(int bit, int npops, Stack * stk, tipId_t pat,
                             int doSing);
static FILE *openOutput(const char *chr);

const char *useMsg =
    "\nUsage: mapmix --admix <fraction> --zero <a.legosim>"
    " --one <b.legosim>\\\n"
    "   <x>=<in_1> <y>=<in_2> ... outgroup=<in_K>\n\n"
    "   where <fraction> is the admixture fraction as estimated by legofit,\n"
    "   <a.legosim> is a file containing the results of a legosim run\n"
    "   with parameters as fitted by legofit, except that one admixture\n"
    "   fraction is set to 0, <b.legosim> is like <a.legosim> except\n"
    "   that this admixture fraction is set to 1, <x> and <y> are arbitrary\n"
    "   labels, and <in_i> are input files in raf format. Labels may not\n"
    "   include the character \":\". Final label must be \"outgroup\".\n"
    "   Writes to a series of files with names like 1.mapmix, X.mapmix,\n"
    "   etc, where \"1\" and \"X\" are the names of chromosomes. If these\n"
    "   files already exist, the program aborts.\n"
    "\n"
    "   If <in_i> file name ends with .gz, input is decompressed using\n"
    "   gunzip.\n";

/// Print usage message and die.
static void usage(void) {
    fputs(useMsg, stderr);
    putc('\n', stderr);
    fprintf(stderr, "   Maximum number of input files: %lu plus outgroup.\n",
            8 * sizeof(tipId_t));
    fputs("\nOptions may include:\n", stderr);
    tellopt("-1 or --singletons", "Use singleton site patterns");
    tellopt("--version", "Print version and exit");
    tellopt("-h or --help", "Print this message");
    exit(1);
}

/// Open a file named <chr>.mapmax for output, or abort if that file
/// already exists or cannot be opened for writing. Return pointer to
/// opened file.
static FILE *openOutput(const char *chr) {
    int status;
    char fname[256];
    status = snprintf(fname, sizeof fname, "%s.mapmix", chr);
    if(status >= sizeof fname) {
        fprintf(stderr,"%s:%d: buffer overflow\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    if( access(fname, F_OK) != -1 ) {
        fprintf(stderr,"%s:%d: ERR: output file %s already exists.\n",
                __FILE__,__LINE__, fname);
        exit(EXIT_FAILURE);
    }
    FILE *fp = fopen(fname, "w");
    if(fp == NULL) {
        fprintf(stderr,"%s:%d: ERR: can't open file %s for writing.\n",
                __FILE__,__LINE__, fname);
        exit(EXIT_FAILURE);
    }
    return fp;
}


/// This stack is local to this file. It provides a bounds-controlled
/// interface to an external array, which is passed as an argument, buff,
/// to Stack_new.
static Stack *Stack_new(int dim, tipId_t buff[dim]) {
    Stack      *self = malloc(sizeof(Stack));
    CHECKMEM(self);
    self->dim = dim;
    self->buff = buff;
    self->nused = 0;
    return self;
}

/// Frees the stack but not the underlying buffer.
static void Stack_free(Stack * stk) {
    free(stk);
}

/// Add an entry to the stack, checking bounds.
static void Stack_push(Stack * self, tipId_t x) {
    if(self->nused == self->dim) {
        fprintf(stderr, "%s:%s:%d ERR: buffer overflow\n",
                __FILE__, __func__, __LINE__);
        exit(EXIT_FAILURE);
    }
    self->buff[self->nused++] = x;
}

/// Call as generatePatterns(0, npops, stk, 0); Recursive function,
/// which generates all legal site patterns and pushes them onto a
/// stack.
static void
generatePatterns(int bit, int npops, Stack * stk, tipId_t pat, int doSing) {
    if(npops >= 8*sizeof(tipId_t)) {
        fprintf(stderr,"%s:%s:%d: %d is too many populations: max is %lu\n",
                __FILE__,__func__,__LINE__,
                npops, 8*sizeof(tipId_t) - 1);
        exit(EXIT_FAILURE);
    }
    const tipId_t unity = 1;
    if(bit == npops) {
        // Recursion stops here. If current pattern is
        // legal, then push it onto the stack. Then return.

        // Exclude patterns with all bits on, or all bits off.
        if(pat == 0 || pat == (unity << npops) - unity)
            return;
        // Exclude singleton patterns unless "doSing" is true.
        if(!doSing && isPow2(pat))
            return;
        Stack_push(stk, pat);
        return;
    }
    tipId_t     on = unity << bit;
    generatePatterns(bit + 1, npops, stk, pat | on, doSing);    // curr bit on
    generatePatterns(bit + 1, npops, stk, pat, doSing); // curr bit off
}

#define BUFFSIZE 1024

int main(int argc, char **argv) {
    int         i, j, status, optndx, done;
    int         doSing = 0;     // nonzero means use singleton site
                                // patterns
    double      admix = -1.0;
    char        errbuff[BUFFSIZE] = { '\0' };
    const char *lgosim[2] = {NULL, NULL};
    StrDblQueue *queue[2] ={NULL, NULL};
    BranchTab   *lgosim_bt[2] = {NULL, NULL};
    LblNdx      lgosim_lndx[2];

    static struct option myopts[] = {
        // {char *name, int has_arg, int *flag, int val}
        {"admix", required_argument, 0, 'a'},
        {"zero", required_argument, 0, 'z'},
        {"one", required_argument, 0, 'o'},
        {"singletons", no_argument, 0, '1'},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'V'},
        {NULL, 0, NULL, 0}
    };

    // command line arguments
    for(;;) {
        char *end;
        i = getopt_long(argc, argv, "a:z:o:hV1", myopts, &optndx);
        if(i == -1)
            break;
        switch (i) {
        case ':':
        case '?':
            usage();
            break;
        case 'a':
            admix = strtod(optarg, &end);
            if(*end != '\0') {
                fprintf(stderr,"Can't parse %s as a double\n", optarg);
                exit(EXIT_FAILURE);
            }
            break;
        case 'z':
            lgosim[0] = optarg;
            queue[0] = StrDblQueue_parseSitePat(lgosim[0]);
            StrDblQueue_normalize(queue[0]);
            status = LblNdx_from_StrDblQueue(lgosim_lndx, queue[0]);
            if(status) {
                fprintf(stderr,"%s:%d: StrDblQueue has a field that's"
                        " too long:\n",
                __FILE__,__LINE__);
                StrDblQueue_print(queue[0], stderr);
                exit(EXIT_FAILURE);
            }
            lgosim_bt[0] = BranchTab_from_StrDblQueue(queue[0], lgosim_lndx);
            break;
        case 'o':
            lgosim[1] = optarg;
            queue[1] = StrDblQueue_parseSitePat(lgosim[1]);
            StrDblQueue_normalize(queue[1]);
            status = LblNdx_from_StrDblQueue(lgosim_lndx+1, queue[1]);
            if(status) {
                fprintf(stderr,"%s:%d: StrDblQueue has a field that's"
                        " too long:\n",
                __FILE__,__LINE__);
                StrDblQueue_print(queue[1], stderr);
                exit(EXIT_FAILURE);
            }
            lgosim_bt[1] = BranchTab_from_StrDblQueue(queue[1], lgosim_lndx+1);
            break;
        case 'V':
            printf("sitepat version %s\n", GIT_VERSION);
            return 0;
        case 'h':
            usage();
            break;
        case '1':
            doSing = 1;
            break;
        default:
            usage();
        }
    }

    if(admix == -1.0) {
        fprintf(stderr,"Missing argument --admix or -a.\n");
        usage();
    }

    if(lgosim[0][0] == '\0') {
        fprintf(stderr, "%s:%d: missing --zero argument\n",
                __FILE__,__LINE__);
        usage();
    }
    if(lgosim[1][0] == '\0') {
        fprintf(stderr, "%s:%d: missing --one argument\n",
                __FILE__,__LINE__);
        usage();
    }
    checkConsistency(lgosim[0], lgosim[1], queue[0], queue[1]);

    queue[0] = StrDblQueue_free(queue[0]);
    queue[1] = StrDblQueue_free(queue[1]);

    if(!LblNdx_equals(lgosim_lndx, lgosim_lndx+1)) {
        fprintf(stderr, "%s:%d: inconsistent LblNdx objects generated"
                " from --zero and --one.\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    // remaining options: input files
    int         n = argc - optind;  // number of input files
    int         m = n-1;            // number excluding outgroup
    if(m < 2)
        usage();

    char       *poplbl[n];
    char       *fname[n];
    RAFReader  *r[n];

    // Number of inputs must be smaller than the number of bits in an object of
    // type tipId_t.
    if(m > 8*sizeof(tipId_t) - 1) {
        fprintf(stderr, "Error: %d .raf files. Max is %lu.\n",
                n, 8*sizeof(tipId_t) - 1);
        usage();
    }

    LblNdx      lndx;
    LblNdx_init(&lndx);

    // Parse remaining arguments, each of which should be of form
    // x=foo, where x is an arbitrary label and foo is the name of an
    // input file. Last label must be "outgroup".
    for(i = 0; i < n; ++i) {
        fname[i] = poplbl[i] = argv[i + optind];
        (void) strsep(fname + i, "=");
        if(fname[i] == NULL
           || poplbl[i] == NULL
           || strlen(poplbl[i]) == 0
           || strlen(fname[i]) == 0 || strchr(poplbl[i], ':') != NULL)
            usage();
        if(i < m)
            LblNdx_addSamples(&lndx, 1, poplbl[i]);
        r[i] = RAFReader_new(fname[i]);
    }
    if(0 != strcmp("outgroup", poplbl[n-1])) {
        fprintf(stderr,"%s:%d: last label is \"%s\""
                " instead of \"outgroup\".\n",
                __FILE__,__LINE__, poplbl[n-1]);
        usage();
    }

    if(!LblNdx_equals(lgosim_lndx, &lndx)) {
        fprintf(stderr, "%s:%d: LblNdx object generated"
                " from --zero and --one\n"
                " is not consistent with .raf files.\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    printf("# mapmix version %s\n", GIT_VERSION);
    printf("# Population labels:\n");
    for(i = 0; i < n; ++i)
        printf("#  %s=%s\n", poplbl[i], fname[i]);

    // make sure labels are all different
    for(i = 1; i < n; ++i) {
        for(j = 0; j < i; ++j) {
            if(0 == strcmp(poplbl[i], poplbl[j])) {
                fprintf(stderr, "ERR: duplicate labels on command line.\n");
                fprintf(stderr, "     duplicated label: %s\n", poplbl[i]);
                exit(EXIT_FAILURE);
            }
        }
    }

    unsigned long npat = (1UL << m) - 2UL;  // number of site patterns
    if(!doSing)
        npat -= m;
    printf("# %s singleton site patterns.\n",
           (doSing ? "Including" : "Excluding"));
    printf("# Number of site patterns: %lu\n", npat);
    tipId_t     pat[npat];

    {
        // Stack is a interface to array "pat".
        Stack      *stk = Stack_new(npat, pat);

        // Put site patterns into array "pat".
        generatePatterns(0, m, stk, 0, doSing);
        Stack_free(stk);
    }

    // Sort site patterns. Major sort is by number of "on" bits,
    // so that doubleton patterns come first, then tripletons, ets.
    // Secondary sort is by order in which labels are listed
    // on the command line.
    qsort(pat, (size_t) npat, sizeof(pat[0]), compare_tipId);
    fflush(stdout);

    // cond_pr[i] is the conditional probability of admixture given
    // site pattern i.
    double cond_pr[npat];
    for(i=0; i < npat; ++i) {
        // Pr of site pattern i given no admixture.
        double p0 = BranchTab_get(lgosim_bt[0], pat[i]);

        // Pr of site pattern i given admixture
        double p1 = BranchTab_get(lgosim_bt[1], pat[i]);

        // Pr of pattern i and no admixture
        double joint_pr0 = p0*(1.0 - admix);

        // Pr of pattern i and admixture
        double joint_pr1 = p1*admix;

        // Bayes's rule
        cond_pr[i] = joint_pr1 / (joint_pr0 + joint_pr1);
    }

    tipId_t union_all_samples = low_bits_on(m);
    char buff[1000], buff2[1000];

    printf("# admix is conditional prob of admixture"
           " given site pattern.\n");
    printf("%15s %15s\n", "SitePat", "admix");
    for(j = 0; j < npat; ++j) {
        if(!doSing && isPow2(pat[j]))
           continue;
        if(pat[j] == union_all_samples)
            continue;
        snprintf(buff2, sizeof(buff2), "%s",
                 patLbl(sizeof(buff), buff, pat[j], &lndx));
        printf("%15s %15.10lf\n", buff2, cond_pr[j]);
    }

    unsigned long nsites = 0, nfixed=0, nbadaa = 0, nbadref=0, nmultalt=0;

    char currchr[128] = { '\0' };
    FILE *ofp = NULL;

    // Iterate through raf files
    RAFReader_clearChromosomes(n, r);
    done=0;
    while( !done ) {
        status = RAFReader_multiNext(n, r);
        if(status==0)
            status = RAFReader_findDaf(n, r);
        switch(status) {
        case 0:
            ++nsites;
            break;
        case EOF:
            done=1;
            continue;
        case REF_MISMATCH:
            ++nsites;
            ++nbadref;
            continue;
        case MULTIPLE_ALT:
            ++nsites;
            ++nmultalt;
            continue;
        case MONOMORPHIC_SITE:
            ++nsites;
            ++nfixed;
            continue;
        case NO_ANCESTRAL_ALLELE:
            ++nsites;
            ++nbadaa;
            continue;
        default:
            // something wrong.
            mystrerror_r(status, errbuff, sizeof errbuff);
            fprintf(stderr,"%s:%d: input error (%s)\n",
                    __FILE__,__LINE__, errbuff);
            exit(EXIT_FAILURE);
        }

        // p and q are frequencies of derived and ancestral alleles
        double      p[m], q[m];
        for(j = 0; j < m; ++j) {
            p[j] = RAFReader_daf(r[j]); // derived allele freq
            q[j] = 1.0 - p[j];
        }

        // Contribution of current snp to each site pattern.  Inner
        // loop considers each bit in current pattern.  If that bit is
        // on, multiply z by the derived allele frequency, p. If
        // that bit is off, multiply by q=1-p. In the end, z is Prod
        // p[j]^bit[j] * q[j]^(1-bit[j]) where bit[j] is the value (0
        // or 1) of the j'th bit in the pattern.
        double pr = 0.0;
        for(i = 0; i < npat; ++i) {
            tipId_t     pattern = pat[i];
            double      z = 1.0;
            for(j = 0; j < m; ++j) {
                if(pattern & 1u) // Is bit j "on" in current pattern?
                    z *= p[j];
                else
                    z *= q[j];
                pattern >>= 1u;
            }
            if(!isfinite(z)) {
                fprintf(stderr, "%s:%d nonfinite z=%lf\n",
                        __FILE__, __LINE__, z);
                fprintf(stderr, "   pattern=%d\n", pat[i]);
                for(j = 0; j < m; ++j)
                    fprintf(stderr, "   %d: p=%lf q=%lf\n", j, p[j], q[j]);
            }
            assert(0 == pattern);

            // z is probability that a random subsample, consisting of
            // one haploid genome from each population would exhibit
            // pattern i at the current nucleotide site. cond_pr[i] is
            // conditional probability of admixture given site pattern
            // i.
            pr += z * cond_pr[i];
        }
        if(0 != strcmp(currchr, r[0]->chr)) {
            // New chromosome. Open new output file
            if(ofp != NULL)
                fclose(ofp);

            status = snprintf(currchr, sizeof currchr, "%s", r[0]->chr);
            if(status >= sizeof currchr) {
                fprintf(stderr,"%s:%d: buffer overflow\n",__FILE__,__LINE__);
                exit(EXIT_FAILURE);
            }
            ofp = openOutput(r[0]->chr);
            fprintf(ofp, "%s\t%s\n", "pos", "admix");
        }
        
        fprintf(ofp, "%lu\t%0.18g\n", r[0]->nucpos, pr);
    }
    if(ofp != NULL)
        fclose(ofp);
    fprintf(stderr, "# Aligned sites                  : %lu\n", nsites);
    if(nbadref)
        fprintf(stderr, "# Disagreements about ref allele : %lu\n", nbadref);
    if(nmultalt)
        fprintf(stderr, "# Sites with multiple alt alleles: %lu\n", nmultalt);
    if(nfixed)
        fprintf(stderr, "# Monomorphic sites              : %lu\n", nfixed);
    if(nbadaa)
        fprintf(stderr, "# Undetermined ancestral allele  : %lu\n", nbadaa);
    fprintf(stderr, "# Sites used                     : %lu\n",
            nsites - nbadaa - nbadref - nmultalt);

    for(i = 0; i < n; ++i)
        RAFReader_free(r[i]);

    fprintf(stderr, "mapmix is finished\n");
    return 0;
}
