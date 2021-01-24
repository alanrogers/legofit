#include "rational.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

long gcd_euclid(long x, long y);

// greatest common divisor: Euclid's method
long gcd_euclid(long x, long y) {
    if(x < 0) x = -x;
    if(y < 0) y = -y;

    if(x==0) return y;
    if(y==0) return x;

    // Euclid's method
    while(x != y) {
        if(x > y)
            x -= y;
        else
            y -= x;
    }

    return x;
}

int main(void) {

    int verbose=0;
    Rational x = Rational_set(105, 140);
    Rational y = Rational_set(18,91);

    long euclid = gcd_euclid(105, 140);
    long lame = gcd(105, 140);

    assert(euclid == lame);
    assert(Rational_eq(x, Rational_set(3,4)));
    assert(Rational_eq(y, Rational_set(36,182)));
    if(verbose) {
        printf("x=");
        Rational_pr(x, stdout);
        printf(", y=");
        Rational_pr(y, stdout);
        putchar('\n');
    }

    Rational z = Rational_add(x, y);
    if(verbose) {
        printf("x+y = ");
        Rational_pr(z, stdout);
        putchar('\n');
    }
    assert(Rational_eq(z, Rational_set(345,364)));

    z = Rational_mul(x, y);
    if(verbose) {
        printf("x*y = ");
        Rational_pr(z, stdout);
        putchar('\n');
    }
    assert(Rational_eq(z, Rational_set(27,182)));

    z = Rational_div(x, y);
    if(verbose) {
        printf("x/y = ");
        Rational_pr(z, stdout);
        putchar('\n');
    }
    assert(Rational_eq(z, Rational_set(91,24)));

    z = Rational_sub(x, y);
    if(verbose) {
        printf("x-y = ");
        Rational_pr(z, stdout);
        putchar('\n');
    }
    assert(Rational_eq(z, Rational_set(201,364)));

    z = Rational_sub(y, x);
    if(verbose) {
        printf("y - x = ");
        Rational_pr(z, stdout);
        putchar('\n');
    }
    assert(Rational_eq(z, Rational_set(-201,364)));

    assert(Rational_eq(Rational_neg(z), Rational_set(201, 364)));
    assert(Rational_eq(Rational_inv(x), Rational_set(4,3)));

    x = Rational_set(LONG_MAX, LONG_MAX);
    assert(Rational_eq(x, Rational_unity));

    // LONG_MIN will overflow
    x = Rational_set(LONG_MIN+1, LONG_MIN+1);
    assert(Rational_eq(x, Rational_unity));

    y = Rational_set(LONG_MAX-1, 1);
    z = Rational_add(x, y);
    assert(z.num == LONG_MAX);
    assert(z.den = 1);

    printf("Rational OK\n");
    
    return 0;
}
