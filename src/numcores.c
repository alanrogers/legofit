#include "misc.h"
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
void usage(void);

void usage(void) {
    fprintf(stderr,"usage: numcores [factor]\n");
    fprintf(stderr,"  where \"factor\" is a number between 0 and 1.\n");
    fprintf(stderr,"  Program prints factor times the number of cores,\n");
    fprintf(stderr,"  rounded to nearest integer.\n");
    fprintf(stderr,"  Default: factor=1.\n");
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
    double factor;
    char *end;
    switch(argc) {
    case 1:
        printf("%d\n", getNumCores());
        break;
    case 2:
        errno=0;
        factor = strtod(argv[1], &end);
        if(errno) {
            fprintf(stderr,"Bad float: %s (%s)\n",
                    argv[1], strerror(errno));
            usage();
        }else if(end == argv[1]) {
            fprintf(stderr,"Bad float: %s\n",
                    argv[1]);
            usage();
        }
        printf("%0.0lf\n", floor(0.5 + factor * getNumCores()));
        break;
    default:
        usage();
        break;
    }
    return 0;
}
