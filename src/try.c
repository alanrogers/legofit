#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>

int main(int argc, char **argv) {

    char buff[100], *leftover;
    double x;

    strcpy(buff, "  1e 4000 ");
    errno=0;
    x = strtod(buff, &leftover);
    if(errno) {
        char err_buff[50];
        strerror_r(errno, err_buff, sizeof(err_buff));
        fprintf(stderr,"%s:%d: Can't parse \"%s\" as float (%s).\n",
                __FILE__,__LINE__, buff, err_buff);
    }
    printf("x=%lf\n", x);

    if(leftover==buff) {
        fprintf(stderr,"%s:%d: Can't parse \"%s\" as float.\n",
                __FILE__,__LINE__, buff);
        exit(1);
    }
    while(isspace(*leftover))
        ++leftover;
    
    if(*leftover != '\0') {
        fprintf(stderr, "%s:%d: Extra chars at end of float \"%s\"."
                " Token=\"%s\".\n",
                __FILE__, __LINE__, leftover, buff);
        exit(1);
    }
	return 0;
}
