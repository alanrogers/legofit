#include <stdio.h>
#include <stdint.h>

int main(void) {
    uint8_t x = 0;

    for(int i=0; i<257; ++i)
        printf("%hhu\n", x++);

    return 0;
}


