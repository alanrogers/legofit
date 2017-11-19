#include <stdio.h>
#include <limits.h>

int main(void) {

    int i = INT_MAX;

    printf("max=%d = %e\n", i, (double) i);

    i += 1;

    printf("max+1=%d\n", i);

    return 0;
}
