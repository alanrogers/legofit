#include <stdio.h>
#include <limits.h>

int main(void) {

    long int i = LONG_MAX;

    printf("max=%ld = %e\n", i, (double) i);

    i += 1;

    printf("max+1=%ld\n", i);

    return 0;
}
