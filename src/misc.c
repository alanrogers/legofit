#include "misc.h"
#include "binary.h"
#include "lblndx.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#ifdef _WIN32
#include <windows.h>
#elif defined(MACOS)
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif
#include <execinfo.h>

#define CALLSTACK_SIZE 128
int         comparePtrs(const void *void_x, const void *void_y);

void dostacktrace(const char *file, int line, FILE * ofp) {
    void       *callstack[CALLSTACK_SIZE];
    int         nsymbols = backtrace(callstack, CALLSTACK_SIZE);

    fprintf(ofp, "backtrace returned %d\n", nsymbols);
    fprintf(ofp, "dostacktrace called from %s:%d:\n", file, line);
    backtrace_symbols_fd(callstack, nsymbols, fileno(ofp));
}

/*
 * Describe an option. For use in "usage" functions.
 */
void tellopt(const char *opt, const char *description) {
    fprintf(stderr, "   %s\n      %s\n", opt, description);
    return;
}

/*
 * An almost platform-independent function that returns the number of
 * CPU cores on the current machine.
 * Source: http://stackoverflow.com/questions/150355/
 * programmatically-find-the-number-of-cores-on-a-machine. If I
 * understand the webpage correctly, this code was written by Dirk-Jan
 * Kroon.
 */
int getNumCores(void) {
#ifdef WIN32
    SYSTEM_INFO sysinfo;

    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#elif defined(MACOS)
    int         nm[2];
    size_t      len = 4;
    uint32_t    count;

    nm[0] = CTL_HW;
    nm[1] = HW_AVAILCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);

    if(count < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &count, &len, NULL, 0);
        if(count < 1) {
            count = 1;
        }
    }
    return count;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

/* eprintf: print error message and exit */
void eprintf(const char *fmt, ...) {
    va_list     args;

    fflush(stdout);

    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);

    if(fmt[0] != '\0' && fmt[strlen(fmt) - 1] == ':')
        fprintf(stderr, " %s", strerror(errno));
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}

/*
 * Uniform perturbation on log scale. Log10 of new value in range
 * (log10(x)-w, log10(x)+w)
 */
double perturb_ratio_w(double x, double w, gsl_rng *rng) {
	assert(x >= 0.0);
	if(x == 0.0)
		return x;

    double lgx = log10(x);
    return pow(10.0, gsl_ran_flat(rng, lgx-w, lgx+w));
}

/*
 * Uniform perturbation on log scale. New value in range (0.1*x, 10*x).
 */
double perturb_ratio(double x, gsl_rng *rng) {
	return perturb_ratio_w(x, 1.0, rng);
}

/*
 * Uniform perturbation within interval. New value in range (lo, hi).
 */
long double perturb_interval(long double x, long double lo, long double hi, gsl_rng *rng) {
    if(!(lo <= x)) {
        fflush(stdout);
        fprintf(stderr,"%s:%s:%d: lo=%Lf must be <= x=%Lf\n",
                __FILE__, __func__, __LINE__, lo, x);
        exit(1);
    }
    if(!(x < hi)) {
        fflush(stdout);
        fprintf(stderr,"%s:%s:%d: x=%Lf must be < hi=%Lf\n",
                __FILE__, __func__, __LINE__, x, hi);
        exit(1);
    }
    long double c_hi = hi - 2 * DBL_EPSILON;
    long double c_lo = lo + 2 * DBL_EPSILON;
    assert(c_lo < c_hi);
	long double rval = gsl_ran_flat(rng, c_lo, c_hi);
	assert(rval > c_lo);
	assert(rval < c_hi);
	return rval;
}

void checkmem( /*@null@ */ void *obj, const char *file, int line) {
    if(obj == NULL)
        eprintf("%s:%d: allocation error", file, line);
}

/** Divide lower triangle of tab by nreps */
void normItab(double tab[3][3], unsigned nreps){
    int i, j;
    for(i = 1; i < 3; ++i) {
        for(j = 0; j < i; ++j) {
            tab[i][j] /= nreps;
        }
    }
}

void printItab(double tab[3][3], FILE *fp){
    int i, j;
    double s=0.0;
    fprintf(fp, "%3s", "");
    for(j = 0; j < 2; ++j)
        fprintf(fp, " %8d", j);
    putc('\n', fp);
#if 1
    for(i = 1; i < 3; ++i) {
        fprintf(fp, "%2d:", i);
        for(j = 0; j < i; ++j) {
            fprintf(fp, " %8.4lf", tab[i][j]);
            s += tab[i][j];
        }
        putc('\n', fp);
    }
#else
    for(i = 1; i < 3; ++i) {
        fprintf(fp, "%2d:", i);
        for(j = 0; j < i; ++j) {
            fprintf(fp, " %0.20le", tab[i][j]);
            s += tab[i][j];
        }
        putc('\n', fp);
    }
#endif
    printf("Sum: %lf\n", s);
}

void printEtab(double x10, double x20, double x21, FILE *fp) {
    int j;
    fprintf(fp, "%3s", "");
    for(j = 0; j < 2; ++j)
        fprintf(fp, " %8d", j);
    putc('\n', fp);
    fprintf(fp, "%2d: %8.4lf\n", 1, x10);
    fprintf(fp, "%2d: %8.4lf %8.4lf\n", 2, x20, x21);
    printf("Sum: %lf\n", x10+x20+x21);
}

/** Print a scatter plot in PicTeX format */
void pictex(double *x, double *y, unsigned n, const char *xlbl,
            const char *ylbl, const char *plotheading,
            const char *fname) {
    double maxval=DBL_MIN, minval=DBL_MAX, unit, range;
    double xaxis = 0.2;
    int i;
    FILE *fp = fopen(fname, "w");
    if(fp == NULL)
        eprintf("%s:%s:%d: can't open \"%s\" for output\n",
                __FILE__, __func__, __LINE__,
                fname);
    for(i=0; i<n; ++i) {
        maxval = fmax(maxval, fmax(x[i], y[i]));
        minval = fmin(minval, fmin(x[i], y[i]));
    }
    range = maxval - minval;
    if(minval < 0.0) {
        maxval += 0.02*range;
        minval -= 0.02*range;
    }else{
        minval = 0.0;
        maxval *= 1.02;
    }
    range = maxval - minval;
    fputs("%-*-latex-*-\n", fp);
    fputs("\\mbox{\\beginpicture\n", fp);
    fputs("\\headingtoplotskip 0pt\n", fp);
    fprintf(fp,"%% Horizontal axis: %0.3lf\\textwidth\n", xaxis);
    unit = xaxis / range;
    fprintf(fp,"\\setcoordinatesystem units"
            " <%0.5lf\\textwidth,%0.5lf\\textwidth>\n",
            unit, unit);
    fprintf(fp,"\\setplotarea x from %0.3lf to %0.3lf,"
            " y from %0.3lf to %0.3lf\n",
            minval, maxval, minval, maxval);
    fputs("\\axis bottom label {Predicted}\n", fp);
    fprintf(fp, "  ticks numbered from %0.3lf to %0.3lf by %0.3lf /\n",
            minval, maxval, range);
    fputs("\\axis left label {\\begin{sideways}Simulated\\end{sideways}}\n", fp);
    fprintf(fp, "  ticks numbered from %0.3lf to %0.3lf by %0.3lf /\n",
            minval, maxval, range);
    fprintf(fp, "\\plotheading{%s}\n", plotheading);
    fprintf(fp, "\\plot %0.3lf %0.3lf %0.3lf %0.3lf /\n",
            minval, minval, maxval, maxval);
    fputs("\\multiput {$\\circ$} at\n", fp);
    fprintf(fp, "%%%6s %6s\n", xlbl, ylbl);
    for(i=0; i<n; ++i)
        fprintf(fp, " %6.4lf %6.4lf\n", x[i], y[i]);
    fputs("/\n", fp);
    fputs("\\endpicture}\n", fp);
    fclose(fp);
}

/** Print a scatter plot in PicTeX format */
void pictex_par(double *x, double *y, unsigned n, const char *xlbl,
            const char *ylbl, const char *plotheading,
            const char *fname) {
    double maxval=DBL_MIN, minval=DBL_MAX, unit, range;
    double max_par = DBL_MIN, min_par=DBL_MAX, par_unit, par_range;
    double xaxis = 0.2;
    int i;
    FILE *fp = fopen(fname, "w");
    if(fp == NULL)
        eprintf("%s:%s:%d: can't open \"%s\" for output\n",
                __FILE__, __func__, __LINE__,
                fname);
    for(i=0; i<n; ++i) {
        maxval = fmax(maxval, x[i]);
        minval = fmin(minval, x[i]);
	max_par = fmax(max_par, y[i]);
	min_par = fmin(min_par, y[i]);
    }
    range = maxval - minval;
    par_range = max_par - min_par;
    if(minval < 0.0) {
        maxval += 0.02*range;
        minval -= 0.02*range;
    }else{
        minval = 0.0;
        maxval *= 1.02;
    }
    if(min_par < 0.0) {
        max_par += 0.02*par_range;
        min_par -= 0.02*par_range;
    }else{
        min_par = 0.0;
        max_par *= 1.02;
    }
    range = maxval - minval;
    par_range = max_par - min_par;
    fputs("%-*-latex-*-\n", fp);
    fputs("\\mbox{\\beginpicture\n", fp);
    fputs("\\headingtoplotskip 0pt\n", fp);
    fprintf(fp,"%% Horizontal axis: %0.3lf\\textwidth\n", xaxis);
    unit = xaxis / range;
    par_unit = xaxis / par_range;
    fprintf(fp,"\\setcoordinatesystem units"
            " <%0.5lf\\textwidth,%0.5lf\\textwidth>\n",
            unit, par_unit);
    fprintf(fp,"\\setplotarea x from %0.3lf to %0.3lf,"
            " y from %0.3lf to %0.3lf\n",
            minval, maxval, min_par, max_par);
    fprintf(fp, "\\axis bottom label %s\n", xlbl);
    fprintf(fp, "  ticks numbered from %0.3lf to %0.3lf by %0.3lf /\n",
            minval, maxval, range);
    fprintf(fp, "\\axis left label {%s}\n", ylbl);
    fprintf(fp, "  ticks numbered from %0.3lf to %0.3lf by %0.3lf /\n",
            min_par, max_par, par_range);
    fprintf(fp, "\\plotheading{%s}\n", plotheading);
    /*
    fprintf(fp, "\\plot %0.3lf %0.3lf %0.3lf %0.3lf /\n",
            minval, min_par, maxval, max_par); */
    fputs("\\multiput {$\\circ$} at\n", fp);
    fprintf(fp, "%%%6s %6s\n", xlbl, ylbl);
    for(i=0; i<n; ++i)
        fprintf(fp, " %6.4lf %6.4lf\n", x[i], y[i]);
    fputs("/\n", fp);
    fputs("\\endpicture}\n", fp);
    fclose(fp);
}

/*
 * Vector v must be sorted in ascending order before this function is
 * called.  Function returns index of first element in sorted array
 * that is >= val.  The function assumes without checking that the
 * input is sorted. If val > vec[len-1], the function returns len.
 */
unsigned Dbl_first_geq(double val, unsigned len, double v[len]) {
    register unsigned lo, mid, hi;

    assert(len > 0);
    lo = 0;
    hi = len - 1;
    if(val > v[hi])
        return len;
    while(lo < hi) {
        mid = lo + (hi - lo) / 2;
        if(mid == lo)
            break;
        if(v[mid] < val)
            lo = mid;
        else
            hi = mid;
    }
    if(v[lo] >= val)
        hi = lo;

    assert(hi >= 0);
    assert(hi < len);
    assert(v[hi] >= val);
    assert(hi == 0 || v[hi - 1] < val);

    return hi;
}

void unitTstResult(const char *facility, const char *result) {
    printf("%-26s %s\n", facility, result);
}

/*
 * In string str, count the number of contiguous chunks of characters
 * belonging to set.
 */
int strCountSetChunks(const char *str, const char *sep) {
    int         nchunks = 0, i;

    while(*str != '\0') {
        i = strcspn(str, sep);  /* skip chars not in set */
        str += i;
        i = strspn(str, sep);   /* skip chars in set */
        if(i > 0) {
            ++nchunks;
            str += i;
        }
    }
    return nchunks;
}

/// duplicate memory block
void       *memdup(const void *p, size_t n) {
    void       *q;

    assert(p != NULL);
    assert(n > 0);

    q = malloc(n);
    checkmem(q, __FILE__, __LINE__);
    memcpy(q, p, n);
    return q;
}

/// Compare pointers to pointers to two tipId_t values.
///
/// @param void_x,void_y pointers to pointers to tipId_t values
/// @returns <0, 0, or >0 depending on whether the first arg is <,
/// ==, or > the second.
int comparePtrs(const void *void_x, const void *void_y) {
    tipId_t * const * x = (tipId_t * const *) void_x;
    tipId_t * const * y = (tipId_t * const *) void_y;

    // Major sort is on the number of samples
    // represented in the site pattern. Patterns with
    // fewer samples come first.
    int diff1bits = num1bits(**x) - num1bits(**y);
    if(diff1bits)
        return diff1bits;

    // Reverse order of bits so that low-order bit
    // is most significant. This ensures that the
    // sort order of samples corresponds to the
    // order in which they were listed in the input
    // data.
    unsigned rx = reverseBits(**x);
    unsigned ry = reverseBits(**y);

    return ry - rx;
}

/// Generate a label for site pattern tid. Label goes into
/// buff. Function returns a pointer to buff;
char       *patLbl(size_t n, char buff[n], tipId_t tid, LblNdx * lblndx) {
    int         maxbits = 40;
    int         bit[maxbits];
    int         i, nbits;
    nbits = getBits(tid, maxbits, bit);
    buff[0] = '\0';
    char        lbl[100];
    for(i = 0; i < nbits; ++i) {
        snprintf(lbl, sizeof(lbl), "%s",
                 LblNdx_lbl(lblndx, (unsigned) bit[i]));
        if(strlen(buff) + strlen(lbl) >= n)
            eprintf("%s:%s:%d: buffer overflow\n", __FILE__, __func__,
                    __LINE__);
        strcat(buff, lbl);
        if(i + 1 < nbits && 1 + strlen(buff) < n)
            strcat(buff, ":");
    }
    return buff;
}

/// On entry, pat is an array of n tipId_t values. On return,
/// ord[0] is the index of the first value, ord[1] is that of the 2nd,
/// and so on.
void orderpat(int n, unsigned ord[n], tipId_t pat[n]) {
    tipId_t  *ptr[n];
    int i;
    for(i=0; i < n; ++i)
        ptr[i] = pat+i;
    qsort(ptr, (size_t) n, sizeof(ptr[0]), comparePtrs);
    for(i=0; i<n; ++i)
        ord[i] = ptr[i]-pat;
}

/// Return Kullback-Leibler divergence of o relative to e.
double KLdiverg(int n, const double o[n], const double e[n]) {
    int i;
    double kl=0.0;
    assert(Dbl_near(1.0, sum_double(n, o)));
    assert(Dbl_near(1.0, sum_double(n, e)));
    for(i = 0; i < n; ++i) {
        if(o[i]==0.0 && e[i]==0.0)
            continue;
        double q = o[i];
        double p = e[i];
        kl += p*log(p/q);
    }
    return kl;
}

/// Sum an array of doubles. Unwrapped loop.
double sum_double(int n, const double x[n]) {
    register double rval;
    register int i, m;

    assert(n >= 0);

    rval = 0.0;
    m = n % 5;
    for(i = 0; i < m; ++i)
        rval += x[i];

    for(i = m; i < n; i += 5)
        rval += x[i] + x[i + 1] + x[i + 2] + x[i + 3] + x[i + 4];

    return rval;
}
