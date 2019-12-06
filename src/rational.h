#ifndef ARR_RATIONAL_H
#define ARR_RATIONAL_H

#include <stdio.h>
typedef struct Rational Rational;

struct Rational {
    long num, den;
};

static const Rational Rational_zero = {.num=0, .den=1};
static const Rational Rational_unity = {.num=1, .den=1};
static const Rational Rational_inf = {.num=1, .den=0};
static const Rational Rational_neginf = {.num=-1, .den=0};
static const Rational Rational_nan = {.num=0, .den=0};

long gcd(long x, long y);
Rational Rational_normalize(Rational x);
Rational Rational_set(long num, long denom);
Rational Rational_add(Rational x, Rational y);
Rational Rational_sub(Rational x, Rational y);
Rational Rational_neg(Rational x);
Rational Rational_inv(Rational x);
Rational Rational_mul(Rational x, Rational y);
Rational Rational_div(Rational x, Rational y);
void Rational_pr(Rational x, FILE *fp);
int Rational_eq(Rational x, Rational y);
long double Rational_ldbl(Rational x);
double Rational_dbl(Rational x);

#endif
