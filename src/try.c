#include <stdio.h>
#include <math.h>
int main(void) {
    __int128_t x = 1234567890987654321;

    printf("%lld%lld\n", (long long) (x>>64), (long long) x);
    return 0;
}
