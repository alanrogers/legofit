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
#include "version.h"
#include "error.h"
#include <assert.h>
#include <ctype.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
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
#ifndef NDEBUG
    void       *callstack[CALLSTACK_SIZE];
    int         nsymbols = backtrace(callstack, CALLSTACK_SIZE);

    fprintf(ofp, "backtrace depth: %d\n", nsymbols);
    fprintf(ofp, "dostacktrace called from %s:%d:\n", file, line);
    backtrace_symbols_fd(callstack, nsymbols, fileno(ofp));
#else
    fprintf(ofp, "Stack trace not available because of NDEBUG option.\n");
#endif
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

/// Open file. On success, return pointer to opened file. On failure
/// abort with an error message. The error message includes the current
/// working directory, whether or not "name" is an absolute pathname.
FILE *efopen(const char *restrict name, const char *restrict mode) {
    FILE *fp = fopen(name, mode);
    if(fp==NULL) {
        char cwd[PATH_MAX], *p;
        p = getcwd(cwd, sizeof(cwd));
        if(NULL == p) {
            fprintf(stderr, "%s:%d: can't open file \"%s.\"\n",
                    __FILE__, __LINE__, name);
        }else{
            fprintf(stderr,"%s:%d: can't open file \"%s.\"\n"
                    "   in directory \"%s\".\n",
                    __FILE__, __LINE__, name, cwd);
        }
        exit(EXIT_FAILURE);
    }
    return fp;
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
    assert(str);
    assert(sep);
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
 * Compare two ints.
 *
 * Function interprets its arguments as pointers to ints.
 *
 * @param void_x,void_y Pointers to the two integers, cast as pointers
 * to voids.
 * @returns -1, 0, or 1 depending on whether the first arg is <,
 * ==, or > the second.
 */
int compareInts(const void *void_x, const void *void_y) {
    const int *x = (const int *) void_x;
    const int *y = (const int *) void_y;

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
    assert(Dbl_near(1.0, sumDbl(n, o)));
    assert(Dbl_near(1.0, sumDbl(n, e)));
    for(i = 0; i < n; ++i) {
        if(o[i]==0.0 && e[i]==0.0)
            continue;
        double q = o[i];
        double p = e[i];
        kl += p*log(p/q);
    }
    return kl;
}

/// Sum an array of doubles.
double sumDbl(int n, const double x[n]) {
    register long double rval = 0.0;

    assert(n >= 0);

    for(register int i = 0; i < n; ++i)
        rval += x[i];

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
#ifdef __clang__
unsigned long strhash(const char *ss) __attribute__((no_sanitize("integer"))){
#else
unsigned long strhash(const char *ss) {
#endif
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
    t = s+1;
    while(*t != '\0')
        ++t;
    while(isspace(*(t-1)))
        --t;
    *t = '\0';
    return s;
}

char* stripInternalWhiteSpace(char* buff) {
  char* str = buff;
  int x = 0;
  for(int i = 0; buff[i]; ++i){
    if(!(buff[i] == ' ')){
      str[x] = str[i];
      ++x;
    }
  }
  str[x] = '\0';
  return str;
}

/**
 * Return tokens separated by 1 or more spaces and/or tabs.
 *
 * On entry, str should be the address of a pointer
 * into a character string. For example:
 *
 *     char buff[100];
 *     char *token, *ptr = buff;
 *
 *     token = nextWhitesepToken(&ptr);
 *
 * On return, token is NULL if no tokens are found. Otherwise, it
 * points to the next token in the string, and *str is either NULL (if
 * the end of the string has been reached) or points to the character
 * immediately following the token.
 */
char *nextWhitesepToken(char **str) {
    char *token;
    do{
        token = strsep(str, " \t\n");
    }while(token!=NULL && *token=='\0');
    return token;
}

/// Separate string s into tokens, separated by any character
/// in string delim. Set ptr[i] equal to the address of the i'th
/// token. Terminate each substring with '\0' within s. Return the
/// number of tokens. Abort if the number of substrings exceeds
/// n.
int tokenize(int dim, char *token[dim], char *s, const char *delim) {
    char *t;
    int n=0;
    while((t=strsep(&s, delim)) != NULL) {
        if(n >= dim) {
            fprintf(stderr,"%s:%d: too many alleles; max is %d\n",
                    __FILE__,__LINE__, dim);
            exit(1);
        }
        token[n] = t;
        ++n;
    }
    return n;
}

/// In string s, replace instances of character a with character b.
void        strReplaceChr(char *s, int a, int b) {
    while(*s != '\0') {
        if(*s == a)
            *s = b;
        ++s;
    }
}

/// Parse token as a double. If strtod sets errno, then leave
/// errno alone. If the token isn't a float or if there are extraneous
/// characters at the end of the token, then set errno=EINVAL.
/// If parseDbl sets errno, it returns 0.0.
double parseDbl(char *token) {
    char *leftover;
    double x;
    errno=0;
    x = strtod(token, &leftover);
    if(errno) {
        // strtod detected a problem
        return 0.0;
    }
    if(leftover==token) {
        // token not interpretable as a float
        errno = EINVAL;
        return 0.0;
    }
    while(isspace(*leftover))
        ++leftover;
    if(*leftover != '\0') {
        // extraneous characters after float
        errno = EINVAL;
        return 0.0;
    }
    return x;
}

/// Truncate string to n characters (plus terminating '\0') by
/// removing characters from the left.
char * strltrunc(char *s, int n) {
    int w = strlen(s);
    if(w <= n)
        return s;

    char *u = s, *v = s + (w-n);
    while(*v)
        *u++ = *v++;
    *u = '\0';
    return s;
}

void hdr(const char *msg) {
    int i, wid, len1, len2, status;
    char version[100], buff[100];
    status = snprintf(version, sizeof(version), "version %s", VERSION);
    if(status >= sizeof(version)) {
        fprintf(stderr,"%s:%d: buffer overflow\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    len1 = strlen(msg);
    len2 = strlen(version);
    wid = 2 + (len1 > len2 ? len1 : len2);
    for(i=0; i < 2 + wid; ++i)
        putchar('#');
    putchar('\n');
    printf("#%s#\n", strcenter(msg, wid, buff, sizeof(buff)));
    printf("#%s#\n", strcenter(version, wid, buff, sizeof(buff)));
    for(i=0; i < 2 + wid; ++i)
        putchar('#');
    putchar('\n');
    putchar('\n');
}

/**
 * Center string "text" in a field of width "width". The centered string
 * is written into the character string "buff", whose size is
 * "buffsize".
 */
char       *strcenter(const char *text, unsigned width,
                      char *buff, size_t buffsize) {
    int         i, j, lpad = 0, rpad = 0, txtwid;

    txtwid = strlen(text);
    if(txtwid >= buffsize) {
        snprintf(buff, buffsize, "%s", text);
        return buff;
    }
    if(width > txtwid)
        lpad = (width - txtwid) / 2;
    rpad = width - txtwid - lpad;
    for(i = 0; i < lpad; ++i)
        buff[i] = ' ';
    for(j = 0; j < txtwid; ++j)
        buff[i + j] = text[j];
    for(j = 0; j < rpad; ++j)
        buff[i + txtwid + j] = ' ';
    buff[lpad + txtwid + rpad] = '\0';
    return buff;
}

/**
 * Remove zeroes from an array of unsigned ints by sliding positive
 * entries to the left. Return number of non-zero entries. This
 * function doesn't re-allocate the array. It just moves zeroes to the
 * end.
 */
int removeZeroes(int dim, unsigned x[dim]) {
    int i, j;
    i=j=0;
    while(i < dim) {
        if(x[j] > 0) {
            // ++ or  0+
            ++i;
            ++j;
        }else if(x[i]==0 && x[j]==0) {
            // 00
            ++i;
        }else if(x[i] > 0 && x[j] == 0) {
            // +0
            assert(i > j);
            x[j] = x[i];
            x[i] = 0;
            ++i;
            ++j;
        }
    }
    while(j<dim && x[j]>0)
        ++j;
    return j;
}

/// Read a line of input into buff.
/// Return EOF or BUFFER_OVERFLOW on failure; 0 on success.
int readline(int dim, char buff[dim], FILE *fp) {
    if(fgets(buff, dim, fp) == NULL)
        return EOF;

    if(NULL == strchr(buff, '\n')) {
        if(feof(fp))
            return 0;
        int c = fgetc(fp);
        if(c == EOF)
            return 0;
        else {
            ungetc(c, fp);
            return BUFFER_OVERFLOW;
        }
    }

    return 0;
}

/// Strip leading pathname components; return a pointer
/// to the filename. There is a standard library function for this,
/// but it can't be used with const strings. This one can.
/// Current code assumes that '/' is the path separator.
/// If pathname ends with '/', this will return an empty string.
const char *mybasename(const char *pathname) {
    const char *p = rindex(pathname, '/');
    if(p == NULL)
        return pathname;
    return p+1;
}

/// Return nonzero if name is legal; zero otherwise. Legal
/// names begin with an alphabetic character and then continue
/// with any number of letters, digits, underscores, colons,
/// or periods".
int legalName(const char *name) {
    const char *legal =
        "abcdefghijklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ" "0123456789" "_:.";
    if(!isalpha(name[0]))
        return 0;
    return strlen(name) == strspn(name, legal);
}

// Copy characters from src to dst until either a '\0' character is
// encountered or m characters have been copied.  Make sure src is
// terminated with '\0'. If dst is too small to hold null-terminated
// string, return BUFFER_OVERFLOW. Otherwise, return 0.
int strnncopy(size_t n, char dst[n], size_t m, const char src[m]) {
    if(n==0)
        return BUFFER_OVERFLOW;
    int i;
    for(i=0; i < m; ++i) {
        if(i >= n) {
            dst[n-1] = '\0';
            return BUFFER_OVERFLOW;
        }
        dst[i] = src[i];
        if(dst[i] == '\0')
            return 0;
    }
    if(i >= n) {
        dst[n-1] = '\0';
        return BUFFER_OVERFLOW;
    }
    dst[i] = '\0';
    return 0;
}

// Count the number of character c in string s. Every string
// is assumed to have exactly one copy of the '\0' character.
int strchrcnt(const char *s, int c) {
    if(c == '\0')
        return 1;
    char *p = strchr(s, c);
    int n=0;
    while(p != NULL) {
        ++n;
        s = p+1;
        p = strchr(s, c);
    }
    return n;
}

// Replace multiple whitespace characters with a single space
// in buff.
void collapse_whitespace(char * buff) {
    char *w, *r;
    int nsp = 0;
    w = r = buff;
    while(*r != '\0') {
        if(isspace(*r)) {
            if(nsp == 0) {
                if(r>w)
                    *w = ' ';
                ++w;
            }
            ++r;
            ++nsp;
        }else{
            nsp = 0;
            if(r>w)
                *w = *r;
            ++r;
            ++w;
        }
    }
    *w = '\0';
}
