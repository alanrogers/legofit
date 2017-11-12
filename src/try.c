#include <stdio.h>

double f(void);

double f(void) {
    return 0.0;
}

int main(void) {

    const void *ptr = f;

    printf("%lf\n", ((double (* const)(void)) ptr)());

    return 0;
}
