/**
 * @file misc.c
 * @author Alan R. Rogers
 * @brief Functions that didn't seem to belong anywhere else.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "misc.h"
#include "binary.h"
#include "lblndx.h"
#include <assert.h>
#include <ctype.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
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

/// Print call stack. This function is not useful when compiler
/// optimizations are turned on, because then the executable
/// doesn't know the names of functions. The stack trace shows
/// addresses instead.
/// Call like this: "dostacktrace(__FILE__, __LINE__)"
void dostacktrace(const char *file, int line, FILE * ofp) {
    void       *callstack[CALLSTACK_SIZE];
    int         nsymbols = backtrace(callstack, CALLSTACK_SIZE);

    fprintf(ofp, "backtrace depth: %d\n", nsymbols);
    fprintf(ofp, "dostacktrace called from %s:%d:\n", file, line);
    backtrace_symbols_fd(callstack, nsymbols, fileno(ofp));
}

/// Describe an option. For use in "usage" functions.
/// @param[in] opt Name of option.
/// @param[in] description Description of option.
void tellopt(const char *opt, const char *description) {
    fprintf(stderr, "   %s\n      %s\n", opt, description);
    return;
}

/**
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

/// eprintf: print error message and exit.
/// This function is dangerous, because unlike ordinary
/// printf, it doesn't detect mismatches in the type and number
/// of arguments.
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

/**
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

/**
 * Uniform perturbation on log scale. New value in range (0.1*x, 10*x).
 */
double perturb_ratio(double x, gsl_rng *rng) {
	return perturb_ratio_w(x, 1.0, rng);
}

/**
 * Uniform perturbation within interval. New value in range (lo, hi).
 */
long double perturb_interval(long double x, long double lo, long double hi,
                             gsl_rng *rng) {
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

/**
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

/**
 * Vector v must be sorted in ascending order before this function is
 * called.  Function returns index of first element in sorted array
 * that is >= val.  The function assumes without checking that the
 * input is sorted. If val > vec[len-1], the function returns len.
 */
long long_first_geq(long val, long *v, long len) {
    register long lo, mid, hi;

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

/**
 * Vector v must be sorted in ascending order before this function is
 * called.  Function returns index of last element in sorted array
 * that is <= val.  The function assumes without checking that the
 * input is sorted. If val < vec[0], the function returns -1.
 */
long long_last_leq(long val, long *v, long len) {
    register long lo, mid, hi;

    assert(len > 0);
    lo = 0;
    hi = len - 1;
    if(val < v[0])
        return -1;
    while(lo < hi) {
        mid = hi - (hi - lo) / 2;
        if(mid == hi)
            break;
        if(v[mid] > val)
            hi = mid;
        else
            lo = mid;
    }
    if(v[hi] <= val)
        lo = hi;

    assert(lo >= 0);
    assert(lo < len);
    assert(v[lo] <= val);
    assert(lo == len - 1 || v[lo + 1] > val);

    return lo;
}

/// Print result of a unit test.
/// @param[in] facility the name of the facility being tested
/// @param[in] result something like "OK" or "FAIL".
void unitTstResult(const char *facility, const char *result) {
    printf("%-26s %s\n", facility, result);
}

/**
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
    assert(p);
    assert(n > 0);
    void       *q;

    q = malloc(n);
    CHECKMEM(q);
    memcpy(q, p, n);
    return q;
}

/**
 * Compare two long ints.
 *
 * Function interprets its arguments as pointers to long ints.
 *
 * @param void_x,void_y Pointers to the two integers, cast as pointers
 * to voids.
 * @returns -1, 0, or 1 depending on whether the first arg is <,
 * ==, or > the second.
 */
int compareLongs(const void *void_x, const void *void_y) {
    const long *x = (const long *) void_x;
    const long *y = (const long *) void_y;

    return (*x < *y) ? -1 : (*x > *y) ? 1 : 0;
}

/**
 * Compare two doubles.
 *
 * Function interprets its arguments as pointers to doubles.
 *
 * @param void_x,void_y Pointers to the two doubles, cast as pointers
 * to voids.
 * @returns -1, 0, or 1 depending on whether the first arg is <,
 * ==, or > the second.
 */
int compareDoubles(const void *void_x, const void *void_y) {
    const double *x = (const double *) void_x;
    const double *y = (const double *) void_y;

    return (*x < *y) ? -1 : (*x > *y) ? 1 : 0;
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

/**
 * Fold x back and forth across the boundaries "lo" and "hi" to obtain a value
 * y such that lo <= y <= hi.
 */
double reflect(double x, double lo, double hi) {
    assert(hi > lo);
    x -= lo;
    hi -= lo;

    double      z = fabs(fmod(x, 2.0 * hi));

    /*    printf("initially z=%lg\n", z); */

    if(z > hi)
        z = 2.0 * hi - z;

    z += lo;

    assert(z >= lo && z <= hi + lo);

    return z;
}

/// Convert NULL-terminated string to lower case
char       *strlowercase(char *s) {
    char *p;

    for(p = s; *p != '\0'; ++p)
        *p = tolower(*p);
    return s;
}

/// Hash a character string
unsigned strhash(const char *ss) {
    unsigned long hashval;
    int c;
    const unsigned char *s = (const unsigned char *) ss;

    // djb2
    hashval = 5381;
    while((c = *s++))
        hashval = ((hashval << 5) + hashval) +  c;

    return hashval;
}

/// Remove character c from string s. Return length of string.
int stripchr(char *s, int c) {
    char *p, *q;
    p = q = s;
    while(*q != '\0') {
        if(*q == c)
            ++q;
        else if(p==q) {
            ++p;
            ++q;
        }else
            *p++ = *q++;
    }
    *p = '\0';
    return p - s;
}

/// Remove white space from beginning and end of a string.
/// Return pointer to stripped string, which lies within the
/// memory buffer passed as input.
char *stripWhiteSpace(char *buff) {
    char *s, *t;
    s = buff;
    while(*s!='\0' && isspace(*s))
        ++s;
    if(*s == '\0')
        return s;
    t = s;
    while(*(t+1) != '\0')
        ++t;
    while(t > s && isspace(*t))
        --t;
    *t = '\0';
    return s;
}

