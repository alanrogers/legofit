#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

void generateSetPartitions(int n_in, int m_in);
void f(int mu, int nu, int sigma);
void b(int mu, int nu, int sigma);

static int n, m;
static int *a = NULL; // array of n+1 ints

void visit(void) {
    for(int i=0; i<=n; ++i) {
        printf("%d", a[i]);
        if(i<n)
            fputs(", ", stdout);
    }
    putchar('\n');
}

/// Generate all ways of partitioning a set of n objects into m
/// non-empty subsets.  This code implements Knuth's answer to
/// exercise 17 of section 7.2.1.5, on pp 764-765 of The Art of
/// Computer Programming, Volume 4A, part 1. On p 417, Knuth
/// attributes the algorithm to Frank Ruskey [Lecture Notes in
/// Comp. Sci. 762 (1993), 205-206].

void generateSetPartitions(int n_in, int m_in) {
    n = n_in;
    m = m_in;
    a = malloc((n+1) * sizeof(a[0]));
    if(a==NULL) {
        fprintf(stderr,"%s:%d: malloc error\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    memset(a, 0, (n+1)*sizeof(a[0]));
    for(int j=1; j <= m; ++j)
        a[n-m+j] = j-1;
    f(m, n, 0);
    free(a);
    a = NULL;
}

// Forward recursion
void f(int mu, int nu, int sigma) {
    if(mu == 2)
        visit();
    else
        f(mu-1, nu-1, (mu+sigma) % 2);
    if( nu == mu + 1) {
        a[mu] = mu - 1;
        visit();
        while(a[nu] > 0) {
            a[nu] -= 1;
            visit();
        }
    }else if(nu > mu+1){
        if( (mu + sigma) % 2 == 1 )
            a[nu-1] = mu-1;
        else
            a[mu] = mu - 1;
        if( (a[nu] + sigma) % 2 == 1)
            b(mu, nu-1, 0);
        else
            f(mu, nu-1, 0);
        while( a[nu] > 0) {
            a[nu] -= 1;
            if( (a[nu] + sigma) % 2 == 1)
                b(mu, nu-1, 0);
            else
                f(mu, nu-1, 0);
        }
    }
}

// Backward recursion
void b(int mu, int nu, int sigma) {
    if( nu == mu+1 ) {
        while( a[nu] < mu-1) {
            visit();
            a[nu] += 1;
        }
        visit();
        a[mu] = 0;
    }else if( nu > mu+1 ) {
        if( (a[nu] + sigma) % 2 == 1)
            f(mu, nu-1, 0);
        else
            b(mu, nu-1, 0);
        while( a[nu] < mu-1 ) {
            a[nu] += 1;
            if( (a[nu] + sigma) % 2 == 1)
                f(mu, nu-1, 0);
            else
                b(mu, nu-1, 0);
        }
        if( (mu + sigma) % 2 == 1)
            a[nu - 1] = 0;
        else
            a[mu] = 0;
    }
    if( mu == 2)
        visit();
    else
        b(mu-1, nu-1, (mu+sigma) % 2);
}
            
int main(void) {

    generateSetPartitions(10, 4);

    return 0;
}
