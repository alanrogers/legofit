/// Partition a positive integer n into a sum of m positive integers.
/// Algorithm H, p 392 of Knuth, Donald E. 2011. The art of computer
/// programming, volume 4A.
#include <stdio.h>
#include <stdlib.h>

void intpart(int n, int m);
void visit(int m, int a[m+1]);

void visit(int m, int a[m+1]) {
    int s=0;
    static int count=0;
    printf("%2d:", ++count);
    for(int i=0; i < m; ++i) {
        printf("%d", a[i]);
        s += a[i];
        if(i < m-1)
            putchar('+');
    }
    printf(" = %d", s);
    putchar('\n');
}

void intpart(int n, int m) {
    if(m < 2) {
        fprintf(stderr,"%s:%d: m (%d) must be > 2\n",
                __FILE__,__LINE__, m);
        exit(EXIT_FAILURE);
    }
    if(m > n) {
        fprintf(stderr,"%s:%d: m (%d) must be <= n (%d)\n",
                __FILE__,__LINE__, m, n);
        exit(EXIT_FAILURE);
    }
    int j;
    int a[m+1];
    a[0] = n - m + 1;
    for(j=1; j < m; ++j)
        a[j] = 1;
    a[m] = -1;

    while(1) {
        visit(m, a);
        if(a[1] < a[0] - 1) {
            a[0] -= 1;
            a[1] += 1;
            continue;
        }

        j = 2;
        int s = a[0] + a[1] - 1;
        while(a[j] >= a[0] - 1) {
            s += a[j];
            j += 1;
        }
        if(j+1>m)
            break;
        int x = a[j] + 1;
        a[j] = x;
        j -= 1;
        while(j>0) {
            a[j] = x;
            s -= x;
            j -= 1;
        }
        a[0] = s;
    }
}

int main(void) {
    int n = 8, m = 8;
    intpart(n, m);
    return 0;
}
