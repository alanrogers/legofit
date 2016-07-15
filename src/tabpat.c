/**
 * @file tabpat.c
 * @brief Tabulate site pattern frequencies from vcf files.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void usage(void);

const char *useMsg =
    "Usage: tabpat <x>=<in1.vcf> <y>=<in2.vcf> ...\n"
    "   where <x> and <y> are arbitrary labels, and <in1> and <in2> are\n"
    "   the names of input files in vcf format. Writes to standard output.\n";

void usage(void) {
    fputs(useMsg, stderr);
    exit(1);
}

int main(int argc, char **argv) {
    int i;
    int n = argc-1; // number of inputs
    FILE *in[n];
    char *poplbl[n];
    char *fname[n];

    if(n == 0)
        usage();

    // Parse arguments, each of which should be of form
    // x=foo, where x is an arbitrary label and foo is the
    // name of an input file.
    for(i=0; i<n; ++i) {
        fname[i] = poplbl[i] = argv[i+1];
        (void) strsep(fname+i, "=");
        if(fname[i] == NULL || strlen(fname[i])==0)
            usage();
        printf("%4s = %s\n", poplbl[i], fname[i]);
        in[i] = fopen(fname[i], "r");
        if(in[i] == NULL) {
            fprintf(stderr,"Couldn't open \"%s\"\n", fname[i]);
            exit(1);
        }
    }

    for(i=0; i<n; ++i)
        fclose(in[i]);
    return 0;
}

