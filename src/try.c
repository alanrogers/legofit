#include <stdio.h>
#include <math.h>
int main(void) {
    double x = -INFINITY;
    if(isinf(x))
        printf("x is infinite\n");
    else
        printf("x is finite\n");
    return 0;
}
