/**
@file rational.c
@page rational
@author Alan R. Rogers
@brief Exact arithmetic with rational numbers.

# `rational`: a library for exact arithmetic with rational numbers

Compile with

    -fsanitize=undefined,integer -fno-sanitize-recover=undefined,integer

@copyright Copyright (c) 2016, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/
#include "rational.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <assert.h>

/// Greatest common divisor: Lame's method
long gcd(long x, long y) {
    if(x < 0) x = -x;
    if(y < 0) y = -y;

    if(x==0) return y;
    if(y==0) return x;

    long r; // remainder
    while(1) {
        if(x > y) {
            r = x % y;
            if(r == 0)
                return y;
            x = r;
        }else{
            r = y % x;
            if(r == 0)
                return x;
            y = r;
        }
    }
    // NOTREACHED
}

/// Remove common factors from numerator and denominator. Make
/// denominator positive.
Rational Rational_normalize(Rational x) {
    if(x.den < 0) {
        x.num = -x.num;
        x.den = -x.den;
    }

    long divisor = gcd(x.num, x.den);
    x.num /= divisor;
    x.den /= divisor;

    return x;
}

/// Set value, normalize.
Rational Rational_set(long num, long den) {
    Rational x;
    x.num = num;
    x.den = den;
    return Rational_normalize(x);
}

/// Add
Rational Rational_add(Rational x, Rational y) {
    Rational z;
    if( x.den == y.den ) {
        z.num = x.num + y.num;
        z.den = x.den;
    }else{
        long l, r, m;
        m = gcd(x.den, y.den);
        l = x.den/m;
        r = y.den/m;
        z.num = (x.num * r) + (y.num * l);
        z.den = x.den * r;
    }
    return Rational_normalize(z);
}

/// Subtract
Rational Rational_sub(Rational x, Rational y) {
    Rational z;
    if( x.den == y.den ) {
        z.num = x.num - y.num;
        z.den = x.den;
    }else{
        z.num = (x.num * y.den) - (y.num * x.den);
        z.den = x.den * y.den;
    }
    return Rational_normalize(z);
}

/// Negate
Rational Rational_neg(Rational x) {
    x.num = -x.num;
    return x;
}

/// Invert
Rational Rational_inv(Rational x) {
    long y = x.num;
    x.num = x.den;
    x.den = y;
    return x;
}

/// Multiply
/// Before multiplying, this function first removes the greatest
/// common divisor of (x.num, y.den) and (x.den, y.num). This reduces
/// the chance of an integer overflow. No comparable reduction is done
/// to the numerator and divisor of x or of y, because I assume they
/// are already in normal form.
Rational Rational_mul(Rational x, Rational y) {
    long divisor = gcd(x.num, y.den);
    x.num /= divisor;
    y.den /= divisor;
    divisor = gcd(x.den, y.num);
    x.den /= divisor;
    y.num /=divisor;
    x.num = x.num * y.num;
    x.den = x.den * y.den;
    return Rational_normalize(x);
}

/// Divide
Rational Rational_div(Rational x, Rational y) {
    x.num = x.num * y.den;
    x.den = x.den * y.num;
    return Rational_normalize(x);
}

/// Print
void Rational_pr(Rational x, FILE *fp) {
    fprintf(fp, "%ld/%ld", x.num, x.den);
}

/// Return 1 if x and y are equal, 0 otherwise.
int Rational_eq(Rational x, Rational y) {
    if(0 == memcmp(&x, &y, sizeof(Rational)))
        return 1;
    return 0;
}

/// Convert to long double
long double Rational_ldbl(Rational x) {
    long double num = x.num;
    long double den = x.den;
    return num/den;
}

/// Convert to double
double Rational_dbl(Rational x) {
    double num = x.num;
    double den = x.den;
    return num/den;
}
