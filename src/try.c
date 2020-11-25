#include <stdio.h>
#include <stdint.h>

int main(void) {
    unsigned u, i, two=2;

    for(i=1; i < 10; ++i) {
        u = 0;
        u = ~u;

        // form 2^i - 1 without forming 2^i
        u >>= 8*sizeof(unsigned) - i - 1;

        printf("%u: u=%u 2^%u-1=%u\n",
               i, u, i, (two<<i) - 1u);
    }

    return 0;
}


