/**
@file pripos.c
@page pripos
@brief Convert prior prob of admixture into posterior prob

**/

#include "branchtab.h"
#include "error.h"
#include "lblndx.h"
#include "misc.h"
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

const char *useMsg =
    "\nUsage: pripos --admix <f> --zero <a.legosim>"
    " --one <b.legosim> <obs.opf>\n\n"
    "   where <f> is the admixture fraction as estimated by legofit,\n"
    "   <a.legosim> is a file containing the results of a legosim run\n"
    "   with parameters as fitted by legofit, except that one admixture\n"
    "   fraction is set to 0, <b.legosim> is like <a.legosim> except\n"
    "   that this admixture fraction is set to 1, and <obs.opf>\n"
    "   is a file of observed site pattern frequencies, as generated\n"
    "   by sitepat.\n\n"
    "Function calculates the conditional probability of admixture for\n"
    "each site pattern and the mean posterior probability of admixture.\n";

/// Print usage message and die.
static void usage(void) {
    fputs(useMsg, stderr);
    putc('\n', stderr);
    fputs("\nOptions may include:\n", stderr);
    tellopt("-1 or --singletons", "Use singleton site patterns");
    tellopt("--version", "Print version and exit");
    tellopt("-h or --help", "Print this message");
    exit(1);
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
    int         i, j, status, optndx;
    int         doSing = 0;     // nonzero means use singleton site
                                // patterns
    double      admix = -1.0;
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
            fprintf(stderr,"%s:%d\n",__FILE__,__LINE__);
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
            fprintf(stderr,"%s:%d\n",__FILE__,__LINE__);
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
            printf("pripos version %s\n", GIT_VERSION);
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

    // There should be one remaining option: the observed.opf
    // input file.
    int n = argc - optind;  // number of remaining arguments
    if(n != 1) {
        fprintf(stderr,"Expecting 1 additional arg: the observed.opf"
                " input. Got %d args instead.\n", n);
        exit(EXIT_FAILURE);
    }

    LblNdx      lndx;
    LblNdx_init(&lndx);

    fprintf(stderr,"%s:%d: parsing %s\n",
            __FILE__,__LINE__, argv[i+optind]);
    StrDblQueue *obs_q = StrDblQueue_parseSitePat(argv[i+optind]);
    fprintf(stderr,"%s:%d\n",__FILE__,__LINE__);
    StrDblQueue_normalize(obs_q);
    status = LblNdx_from_StrDblQueue(&lndx, obs_q);
    if(status) {
        fprintf(stderr,"%s:%d: StrDblQueue has a field that's too long:\n",
                __FILE__,__LINE__);
        StrDblQueue_print(obs_q, stderr);
        exit(EXIT_FAILURE);
    }
    BranchTab *obs_bt = BranchTab_from_StrDblQueue(obs_q, &lndx);
    if(obs_bt == NULL) {
        fprintf(stderr,"%s:%d: can't make obs_bt\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    if(!LblNdx_equals(lgosim_lndx, &lndx)) {
        fprintf(stderr, "%s:%d: LblNdx object generated"
                " from --zero and --one\n"
                " is not consistent with observed data file.\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    printf("# pripos version %s\n", GIT_VERSION);

    unsigned long npat = (1UL << lndx.n) - 2UL;  // number of site patterns
    if(!doSing)
        npat -= lndx.n;
    printf("# %s singleton site patterns.\n",
           (doSing ? "Including" : "Excluding"));
    printf("# Number of site patterns: %lu\n", npat);
    tipId_t     pat[npat];

    {
        // Stack is a interface to array "pat".
        Stack      *stk = Stack_new(npat, pat);

        // Put site patterns into array "pat".
        generatePatterns(0, lndx.n, stk, 0, doSing);
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

    tipId_t union_all_samples = low_bits_on(lndx.n);
    char buff[1000], buff2[1000];
    double posterior_pr = 0.0;

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
        double patfrq = BranchTab_get(obs_bt, pat[j]);
        double joint = patfrq * cond_pr[j];
        posterior_pr += joint;
        printf("%15s %15.10lf\n", buff2, cond_pr[j]);
    }
    fflush(stdout);


    printf("Prior probability of admixture: %lf\n", admix);
    printf("Mean posterior probability of admixture: %lf\n", posterior_pr);

    return 0;
}
