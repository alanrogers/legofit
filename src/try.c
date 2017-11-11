#include <stdio.h>
#include <string.h>
#include <errno.h>

typedef struct Foo {
    const char *name;
    const void *address;
    int type;
    void *context;
} Foo;

int main(int argc, char **argv) {
    double x=2.0, y=3.0;

    Foo vars[] = {{"x", &x}, {"y", &y}};

    printf("%d: name=%s address=%p type=%d\n",
           0, vars[0].name, vars[0].address, vars[0].type);
    printf("%d: name=%s address=%p\n", 1, vars[1].name, vars[1].address);

    return 0;
}
