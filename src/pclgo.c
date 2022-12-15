/**
@file pclgo.c
@page pclgo
@author Alan R. Rogers
@brief Re-express free variables in terms of principal components.

# pclgo: re-express free variables in terms of principal components.

First run `legofit` on a series of bootstrap or simulation replicates
to produce a series of .legofiles. All these `legofit` runs should use
the same .lgo file. `pclgo` reads this .lgo file and also all the
.legofit files. It uses the .legofit files to do a principal
components analysis of the free variables, and then rewrites a portion
of the .lgo file, re-expressing all the free variables as linear
functions of the principle components. In the rewrite, all the free
variables are principal component scores and are therefore
uncorrelated. The output of pclgo is not a complete .lgo file. It must
be integrated into a new .lgo file.

You can do this with a text editor. But here is a trick that uses the
unix shell:

    (grep ^# a.lgo; pclgo a.lgo a2.legofit a2boot*.legofit;\
     grep -v ^# a.lgo | egrep -v "\<free\>") > b.lgo

This assumes that you have used `a.lgo` to analyze a real data set and
several bootstrap replicates. The legofit output files are named
`a2.legofit` (for the real data), and `a2boot*.legofit` (for the
bootstrap replicates). The parentheses group the pipeline within it so
that we can redirect the output with a single `>`. Within the
parentheses, the first command pulls the comments out of a.lgo and
writes then to standard output. They will end up at the top of
`b.lgo`. Next, the `pclgo` command writes lines in `.lgo` format,
which define variables called "pc1", "pc2", etc., which refer to
principal components. Following these, `pclgo` writes lines that
re-express all the free variables of `a.lgo` as linear functions of
the principal components. Following the `pclgo` command is a `grep -v`
command that excludes comments. Its output is piped to an `egrep`
command that excludes the free variables of `a.lgo`. All of this gets
written into `b.lgo`. The result is a .lgo file that re-expresses all
the free variables of `a.lgo` as functions of principal components.

By default, `pclgo` uses all principal components and therefore does
not reduce the dimension of the search space. To reduce dimension, try
`pclgo --tol 0.001`. This would exclude principal components that
account for less than a fraction 0.001 of the variance.

@copyright Copyright (c) 2018, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "strdblqueue.h"
#include "misc.h"
#include "parstore.h"
#include "gptree.h"
#include "network.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

void usage(void);
void parseError(char *p, char *buff, const char *file, int line);
char *lastNonwhite(char *s);

const char *usageMsg =
     "Usage: pclgo [options] <file.lgo> <rep1.legofit> <rep2.legofit> ...\n"
     "  where the .legofit files were generated by legofit, using data sets\n"
     "  corresponding to different bootstrap or simulation replicates. All\n"
     "  these legofit runs should use the same .lgo file, whose name should\n"
     "  be listed first on the pclgo command line.\n"
     "Options:\n"
     "  -p <d>          : specify digits of precision\n"
     "  --tol <f>       : drop dimensions that account for a fraction\n"
     "                    of variance smaller than <f>.\n"
     "  -h or --help    : print this message\n";

void usage(void) {
    fputs(usageMsg, stderr);
    exit(EXIT_FAILURE);
}

void parseError(char *p, char *buff, const char *file, int line) {
    fprintf(stderr,"%s:%d: parse error\n",file, line);
    fprintf(stderr, "  %s\n", buff);
    fprintf(stderr, "  %*s^\nError near here\n", (int)(p-buff), "");
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv){
    int i, j;
    int ncases=0;
    int digits = 16;
    double minFracVar = 0.0;
    const char *letters = "abcdefghijklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    // Command line arguments specify file names
    if(argc < 3)
        usage();

    // first pass through arguments counts file names
    for(i=1; i<argc; ++i) {
        if(0==strcmp("-p", argv[i])) {
            if(++i >= argc)
                usage();
            digits = strtol(argv[i], NULL, 10);
        }else if(0==strcmp("--tol", argv[i])) {
            if(++i >= argc)
                usage();
            minFracVar = strtod(argv[i], NULL);
            if(minFracVar < 0.0 || minFracVar > 1.0) {
                fprintf(stderr,"Argument of --tol not in [0,1]\n");
                exit(EXIT_FAILURE);
            }
        }else if(0==strcmp("-h", argv[i]) || 0==strcmp("--help", argv[i]))
            usage();
        else if(argv[i][0] == '-') {
            fprintf(stderr,"Unknown argument: %s\n", argv[i]);
            usage();
        }else
            ++ncases;
    }
    if(ncases < 2)
        usage();

    fprintf(stderr,"pclgo: %d input files\n", ncases);

    // 2nd pass builds array of bootstrap filenames
    ncases -= 1;   // counts number of bootstrap files
    const char *lgofname = NULL;
    const char *repfname[ncases];
    int gotLgoFile=0;
    for(i=1, j=0; i<argc; ++i) {
        if(0==strcmp("-p", argv[i])) {
            ++i;
            continue;
        }else if(0==strcmp("--tol", argv[i])) {
            ++i;
            continue;
        }else if(argv[i][0] == '-')
            continue;
        else if(!gotLgoFile) {
            lgofname = argv[i];
            gotLgoFile=1;
        }else
            repfname[j++] = argv[i];
    }

    if(lgofname == NULL)
        usage();

    // Read bootstrap files into an array of FIFO queues
    StrDblQueue *queue[ncases];
    for(i=0; i < ncases; ++i) {
        queue[i] = StrDblQueue_parseLegofit(repfname[i]);
        if( queue[i] == NULL ) {
            fprintf(stderr,"%s:%d: can't parse legofit file \"%s\".\n",
                    __FILE__,__LINE__, repfname[i]);
            exit(EXIT_FAILURE);
        }
        if(i>0) {
            if(StrDblQueue_compare(queue[0], queue[i])) {
                fprintf(stderr, "%s:%d: inconsistent parameters in"
                        " files %s and %s\n", __FILE__,__LINE__,
                        repfname[0], repfname[i]);
                exit(EXIT_FAILURE);
            }
        }
    }

    // Use the queues to populate an array of names of free parameters
    // and a matrix of parameter values. Rows are replicates; columns
    // are parameters.
    int npar = StrDblQueue_length(queue[0]);
    char *parname[npar];
    int partype[npar];
    memset(partype, 0, sizeof(partype));
    gsl_matrix *datmat = gsl_matrix_alloc(ncases, npar);
    if(datmat==NULL)
        DIE("bad gsl_matrix_alloc");
    for(i=0; i < ncases; ++i) {
        for(j=0; j < npar; ++j) {
            StrDbl strdbl;
            queue[i] = StrDblQueue_pop(queue[i], &strdbl);
            gsl_matrix_set(datmat, i, j, strdbl.val);
            if(i==0)
                parname[j] = strdup(strdbl.str);
        }
        assert(queue[i] == NULL); // check that queues are freed
    }
    
    // Parse .lgo file in order to check agreement between free
    // variables in the .lgo file and those in the .legofit files.
    Network_init(STOCHASTIC);
    {
        Bounds bnd = {
            .lo_twoN = 1.0,
            .hi_twoN = 1e8,
            .lo_t = 0,
            .hi_t = 1e8
        };
        GPTree *gptree = GPTree_new(lgofname, bnd);
        if(gptree==NULL)
            DIE("bad GPTree_new");
        if(npar != GPTree_nFree(gptree)) {
            fprintf(stderr, "%s:%d: mismatch in count of free parameters.\n",
                    __FILE__,__LINE__);
            fprintf(stderr, "   %s has %d; %s has %d\n",
                    lgofname, GPTree_nFree(gptree),
                    repfname[0], npar);
            GPTree_printParStoreFree(gptree, stderr);
            exit(EXIT_FAILURE);
        }

        for(i=0; i<npar; ++i) {
            const char *lgoParName = GPTree_getNameFree(gptree,i);
            int c = strcmp(parname[i], lgoParName);
            if(c) {
                fprintf(stderr,"%s:%d: mismatch in free parameter names.\n",
                        __FILE__,__LINE__);
                fprintf(stderr," Param %d is %s in file %s but %s in %s.\n",
                        i, lgoParName, lgofname, parname[i], repfname[0]);
                exit(EXIT_FAILURE);
            }
        }
        GPTree_free(gptree);
    }

    double x, y, mean[npar], sd[npar];
    memset(mean, 0, npar*sizeof(mean[0]));
    memset(sd, 0, npar*sizeof(sd[0]));

    // Calculate means
    for(i=0; i < ncases; ++i) {
        for(j=0; j<npar; ++j)
            mean[j] += gsl_matrix_get(datmat, i, j);
    }
    for(j=0; j<npar; ++j) {
        mean[j] /= ncases;
        if( !isfinite(mean[j]) ) {
            fprintf(stderr,"%s:%d: mean[%d] is not finite\n",
                    __FILE__,__LINE__, j);
            exit(EXIT_FAILURE);
        }
    }

    // Calculate standard deviations
    for(i=0; i < ncases; ++i) {
        for(j=0; j<npar; ++j) {
            x = gsl_matrix_get(datmat, i, j);
            x -= mean[j];
            sd[j] += x*x;
        }
    }
    for(j=0; j<npar; ++j) {
        sd[j] = sqrt(sd[j]/(ncases-1));
        if( !isfinite(sd[j]) ) {
            fprintf(stderr,"%s:%d: sd[%d] is not finite\n",
                    __FILE__,__LINE__, j);
            exit(EXIT_FAILURE);
        }
    }

    for(j=0; j < npar; ++j)
        fprintf(stderr,"%10s: mean=%lg, sd=%lg\n", parname[j], mean[j], sd[j]);

    // Rescale data matrix
    for(i=0; i < ncases; ++i) {
        for(j=0; j<npar; ++j) {
            x = gsl_matrix_get(datmat, i, j);
            x = (x - mean[j])/sd[j];
            gsl_matrix_set(datmat, i, j, x);
        }
    }

    // Singular value decomposition.
    int status;
    gsl_matrix *V =
        gsl_matrix_alloc(npar, npar); // eigenvectors
    if(V==NULL)
        DIE("bad gsl_matrix_alloc");
    // s holds singular values
    gsl_vector *s = gsl_vector_alloc(npar);
    if(s == NULL)
        DIE("bad gsl_vector_alloc");

    // Matrix U starts as a copy of datmat and will
    // end up as the scores matrix.
    gsl_matrix *U = gsl_matrix_alloc(ncases, npar);
    if(U==NULL)
        DIE("bad gsl_matrix_alloc");
    status = gsl_matrix_memcpy(U, datmat);
    if(status) {
        fprintf(stderr,"%s:%d: gsl_matrix_memcpy returned %d\n",
                __FILE__,__LINE__, status);
        fprintf(stderr,"    error: %s\n",
                gsl_strerror(status));
        exit(EXIT_FAILURE);
    }

    //#define SVD_MOD
#ifdef SVD_JACOBI
    const char *svdname = "gsl_linalg_SV_decomp_jacobi";
    status = gsl_linalg_SV_decomp_jacobi(U, V, s);
#elif defined( SVD_MOD )
    const char *svdname = "gsl_linalg_SV_decomp_mod";
    {
        gsl_matrix *X = gsl_matrix_alloc(npar, npar);
        if(U==NULL)
            DIE("bad gsl_matrix_alloc");
        gsl_vector *work = gsl_vector_alloc(npar);
        if(work == NULL)
            DIE("bad gsl_vector_alloc");
        status = gsl_linalg_SV_decomp_mod(U, X, V, s, work);
        gsl_matrix_free(X);
        gsl_vector_free(work);
    }
#else
    const char *svdname = "gsl_linalg_SV_decomp";
    {
        gsl_vector *work = gsl_vector_alloc(npar);
        if(work == NULL)
            DIE("bad gsl_vector_alloc");
        status = gsl_linalg_SV_decomp(U, V, s, work);
        gsl_vector_free(work);
    }
#endif
    if(status) {
        fprintf(stderr,"%s:%d: %s returned %d\n",
                __FILE__,__LINE__, svdname, status);
        fprintf(stderr,"    error: %s\n",
                gsl_strerror(status));
        exit(EXIT_FAILURE);
    }

    printf("# PCA calculated with %s\n", svdname);

    // eigenvalues, re-expressed as a fraction of the trace
    double fracvar[npar];  // eigenvalues/trace
    double trace=0.0;   // sum of eigenvalues
    for(i=0; i < npar; ++i) {
        x = gsl_vector_get(s, i);
        if( !isfinite(x) ) {
            fprintf(stderr,"%s:%d: %lf = singular value %d is not finite\n",
                    __FILE__,__LINE__, x, i);
            exit(EXIT_FAILURE);
        }
        fracvar[i] = (x*x)/(ncases-1);
        trace += fracvar[i];
    }

    if( !(trace > 0.0) ) {
        fprintf(stderr,"%s:%d: trace (=%lg) is not positive.\n",
                __FILE__,__LINE__, trace);
            exit(EXIT_FAILURE);
    }
    int nkeep=0;
    for(i=0; i < npar; ++i) {
        fracvar[i] /= trace;
        if(fracvar[i] >= minFracVar)
            ++nkeep;
    }
    printf("# Fraction of variance:\n");
    for(i=0; i<npar; i += 8) {

        putchar('#');
        for(j=0; j<8 && i+j<npar; ++j) {
            char lbl[10];
            status = snprintf(lbl, sizeof(lbl), "pc%d", i+j+1);
            if(status >= sizeof lbl) {
                fprintf(stderr,"%s:%d: buffer overflow\n",__FILE__,__LINE__);
                exit(EXIT_FAILURE);
            }
            printf(" %8s", lbl);
        }
        fputs("\n#", stdout);
        for(j=0; j<8 && i+j<npar; ++j) {
            printf(" %8.5lf", fracvar[i+j]);
        }
        putchar('\n');
    }

    // Calculate score of each PC for each case, and use these
    // to set pcscale, which records the max absolute value of each pc
    // score. It is used to set the range of pc values in the .lgo
    // file.
    double pcscale[npar];
    for(j=0; j < npar; ++j) {
        y = gsl_vector_get(s, j);
        double m = -INFINITY;
        for(i=0; i < ncases; ++i) {
            x = gsl_matrix_get(U, i, j);
            m = fmax(m, fabs(x*y)); // x*y is score for pc j, case i.
        }
        // round up to a power of 2
        pcscale[j] = pow(2.0, ceil(log(m)/M_LOG2E));
    }

    // Rows correspond to original variables, columns to eigenvectors.
    // coef[i][j] is the coefficient of the i'th variable on the j'th
    // eigenvector, rexpressed in units of the original variable.
    double coef[npar][npar];
    for(i=0; i < npar; ++i) {
        for(j=0; j < npar; ++j) {
            x = gsl_matrix_get(V, i, j);
            coef[i][j] = x*sd[i];
        }
    }

    // Parse .lgo file to get the type (twoN, time, etc) of
    // each free variable.
    char buff[500];
    FILE *fp = fopen(lgofname,"r");
    if(fp == NULL) {
        fprintf(stderr,"%s:%d: can't read file \"%s\"\n",
                __FILE__,__LINE__, lgofname);
        exit(EXIT_FAILURE);
    }
    while(1) {
        if(fgets(buff, sizeof buff, fp) == NULL) {
            break;
        }
        if(strchr(buff, '\n') == NULL && !feof(stdin)) {
            fprintf(stderr, "%s:%d: Buffer overflow. size=%zu\n",
                    __FILE__, __LINE__, sizeof(buff));
            exit(EXIT_FAILURE);
        }
        char *p = buff, *q;
        assert(p);
        while(isspace(*p) && *p != '\0')
            ++p;
        if(*p == '#' || *p == '\0')
            continue;
        char tok0[50], tok1[50], varname[50];

        // token 0
        int len = strspn(p, letters);
        if(len == 0)
            parseError(p, buff, __FILE__, __LINE__);
        if(len >= sizeof(tok0))
            DIE("buffer overflow");
        snprintf(tok0, len+1, "%s", p);

        // token 1
        p += len;
        while(isspace(*p) && *p != '\0')
            ++p;
        len = strspn(p, letters);
        if(len == 0)
            parseError(p, buff, __FILE__, __LINE__);
        if(len >= sizeof(tok1))
            DIE("buffer overflow");
        snprintf(tok1, len+1, "%s", p);

        // variable name
        q = strchr(p, '=');
        if(q == NULL)
            continue;
        while(isspace(*(q-1)))
            --q;
        assert(q>p);
        *q = '\0';
        for(p=q-1; p>buff && !isspace(*(p-1)); --p)
            ;
        snprintf(varname, sizeof(varname), "%s", p);

        // We're only interested in free variables
        if(strcmp(tok1, "free") != 0)
            continue;
        for(i=0; i<npar; ++i) {
            // find index of variable
            if(strcmp(varname, parname[i]) == 0)
                break;
        }
        assert(i < npar);
        if(strcmp(tok0, "time") == 0)
            partype[i] = TIME;
        else if(strcmp(tok0, "twoN") == 0)
            partype[i] = TWON;
        else if(strcmp(tok0, "mixFrac") == 0)
            partype[i] = MIXFRAC;
        else if(strcmp(tok0, "param") == 0)
            partype[i] = ARBITRARY;
        else {
            fprintf(stderr,"%s:%d: unknown parameter type: %s\n",
                    __FILE__,__LINE__,tok0);
            exit(EXIT_FAILURE);
        }
    }

    // Output. This is where we put the mean back in.
    for(i=0; i < nkeep; ++i)
        printf("param free [%6lg, %6lg] pc%d = 0\n",
               -pcscale[i], pcscale[i], i+1);
    for(i=0; i < npar; ++i) {
        int pos=0;
        char operator;
        switch(partype[i]) {
        case TIME:
            pos += printf("time constrained");
            break;
        case TWON:
            pos += printf("twoN constrained");
            break;
        case MIXFRAC:
            pos += printf("mixFrac constrained");
            break;
        case ARBITRARY:
            pos += printf("param constrained");
            break;
        default:
            fprintf(stderr,"%s:%d: bad type (%d) for param %s\n",
                    __FILE__,__LINE__, partype[i], parname[i]);
        }
        pos += printf(" %s = %0.*lg", parname[i], digits, mean[i]);
        for(j=0; j < nkeep; ++j) {
            if(coef[i][j] >= 0.0) {
                operator = '+';
                x = coef[i][j];
            }else{
                operator = '-';
                x = -coef[i][j];
            }
            if(pos > 60-digits) {
                printf(" %c\n", operator);
                pos = printf("    %0.*lg*pc%d", digits, x, j+1);
            }else
                pos += printf(" %c %0.*lg*pc%d",
                              operator, digits, x, j+1);
        }
        putchar('\n');
    }

    gsl_matrix_free(datmat);
    gsl_matrix_free(U);
    gsl_matrix_free(V);
    gsl_vector_free(s);
    for(i=0; i < npar; ++i)
        free(parname[i]);
    fclose(fp);
    return 0;
 }
