#include <stdio.h>
int main(void) {

    double x;
    int status;
    while(1) {
        status = scanf("%lf", &x);
        if(status != 1)
            break;
        printf("%g\n", x);
    }

    return 0;
}
