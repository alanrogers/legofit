#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(void) {

    long h;
    char token[100], *end;
    strcpy(token, "-1");
    
    h = strtol(token, &end, 10);
    if(end==token || h<0) // token isn't a nonnegative integer
        printf("%s is Not a nonnegative integer\n", token);
    else           // token is a nonnegative integer
        printf("%s IS a nonnegative integer: value=%ld\n", token, h);

    strcpy(token, " -eI");
    h = strtol(token, &end, 10);
    if(end==token || h<0) // token isn't a nonnegative integer
        printf("%s is Not a nonnegative integer\n", token);
    else           // token is a nonnegative integer
        printf("%s IS a nonnegative integer: value=%ld\n", token, h);

    strcpy(token, " -1");
    h = strtol(token, &end, 10);
    if(end==token || h<0) // token isn't a nonnegative integer
        printf("%s is Not a nonnegative integer\n", token);
    else           // token is a nonnegative integer
        printf("%s IS a nonnegative integer: value=%ld\n", token, h);

    strcpy(token, " 123 ");
    h = strtol(token, &end, 10);
    if(end==token || h<0) // token isn't a nonnegative integer
        printf("%s is Not a nonnegative integer\n", token);
    else           // token is a nonnegative integer
        printf("%s IS a nonnegative integer: value=%ld\n", token, h);
}
