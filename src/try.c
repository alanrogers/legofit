#include <stdio.h>
#include <string.h>
#include <stdlib.h>


int main(int argc, char **argv) {

    int i, n, dim=10;
    char s[100];
    char *tokens[dim];

    strcpy(s, "aaa,bbb,c");
    n = tokenize(dim, tokens, s, ",");
    printf("%d tokens:", n);
    for(i=0; i<n; ++i)
        printf(" `%s'", tokens[i]);
    putchar('\n');

    strcpy(s, "aaa,bbb,c,,");
    n = tokenize(dim, tokens, s, ",");
    printf("%d tokens:", n);
    for(i=0; i<n; ++i)
        printf(" `%s'", tokens[i]);
    putchar('\n');

	return 0;
}
