#include <stdio.h>
#include "binary.h"
int main(void) {
    const uint32_t unity = 1u;
    uint32_t pat;

    pat = ~0u;
    printBits(sizeof pat, &pat, stdout);
    return 0;
}
