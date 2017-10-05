
/**
@file daf.c
@page daf
@brief Calculate derived allele frequency, daf.

Input file should consist of tab-separated columns:

1. chromosome
2. position
3. reference allele
4. alternate alleles
5. ancestral allele
6. genotype in format "0/1" or "0|1", where 0 represents
a copy of the reference allele and 1 a copy of the derived
allele.
7. etc for as many columns as there are genotypes.

With 1000-genomes data, this input can be generated from a vcf or
bcf file as follows:

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA[\t%GT]\n' fname.bcf

Output is in 5 columns, separated by whitespace:

1. chromosome
2. position of the nucleotide
3. aa, the ancestral allele
4. da, the derived allele
5. daf, derived allele frequency

Sites are rejected unless they have a single ref or ancestral
allele. Missing values are allowed for the alt allele. At the end of
the job a summary of rejected sites is written to stderr.

@copyright Copyright (c) 2016, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int main(int argc, char **argv) {

    if(argc != 1) {
        fprintf(stderr, "Usage: daf\n");
        fprintf(stderr, "       Reads standard input;"
                " writes to standard output.\n");
        exit(EXIT_FAILURE);
    }

    const int   buffsize = 4096;
    char        buff[buffsize];

    // Keep track of the number of sites at which the number
    // of reference, alternate, or ancestral alleles differs from 0
    int         zeroref = 0, zeroalt = 0, zeroaa = 0, zerogtype = 0;
    int         missref = 0, missaa = 0;
    int         multref = 0, multalt = 0, multaa = 0;
    int         nbad = 0, ngood = 0;
    int         ok;             // is current line acceptable
    long unsigned lastnucpos=0, nucpos;
    char lastchr[100] = {'\0'};

    printf("#%3s %10s %2s %2s %20s\n", "chr", "pos", "aa", "da", "daf");
    while(1) {
        if(NULL == fgets(buff, buffsize, stdin)) {
            break;
        }
        if(NULL == strchr(buff, '\n') && !feof(stdin)) {
            fprintf(stderr, "%s:%d: Buffer overflow. size=%zu\n",
                    __FILE__, __LINE__, sizeof(buff));
            exit(EXIT_FAILURE);
        }
        char       *chr, *pos, *ref, *alt, *aa, *gtype, *next = buff;

        chr = strsep(&next, "\t");  // field 0
        pos = strsep(&next, "\t");  // field 1
        ref = strsep(&next, "\t");  // field 2
        alt = strsep(&next, "\t");  // field 3
        aa = strsep(&next, "\t");   // field 4

        if(aa == NULL) {
            fprintf(stderr, "%s:%d: Bad input line\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        // lowercase alleles
        strlowercase(ref);
        strlowercase(alt);
        strlowercase(aa);

        // strip extraneous characters
        (void) stripchr(chr, ' ');
        (void) stripchr(pos, ' ');
        int         nref = stripchr(ref, ' ');  // nref: num ref alleles
        int         nalt = stripchr(alt, ' ');  // nalt: num alt alleles
        (void) stripchr(aa, ' ');
        int         naa = stripchr(aa, '|');    // naa: num ancestral alleles
        nucpos = strtoul(pos, NULL, 10);

        // Check sort of chromosomes
        if(*lastchr) {
            int diff = strcmp(lastchr, chr);
            if(diff > 0) {
                // bad sort
                fprintf(stderr,"%s:%d: unsorted chromosomes\n",
                        __FILE__,__LINE__);
                fprintf(stderr,"    %s > %s\n", lastchr, chr);
            }else if(diff < 0) {
                // new chromosome
                int status = snprintf(lastchr, sizeof lastchr, "%s", chr);
                if(status >= sizeof lastchr) {
                    fprintf(stderr,"%s:%d: buffer overflow\n",
                            __FILE__,__LINE__);
                    exit(1);
                }
                lastnucpos = 0;
            }
            assert(diff==0);
        }

        // Check sort of nucleotice positions
        if(lastnucpos) {
            if(lastnucpos >= nucpos) {
                fprintf(stderr,"%s:%d: Duplicate: chr=%s pos=%lu\n",
                        __FILE__,__LINE__, chr, nucpos);
                fprintf(stderr,"%s:%d: Previous : chr=%s pos=%lu\n",
                        __FILE__,__LINE__, lastchr, lastnucpos);
                exit(1);
            }
        }
        lastnucpos = nucpos;

        // Skip sites at which the number of reference, alternate, or ancestral
        // alleles differs from 1.
        ok = 1;
        if(nref == 0) {
            ++zeroref;
            ok = 0;
        }
        if(nalt == 0) {
            ++zeroalt;
            ok = 0;
        }
        if(naa == 0) {
            ++zeroaa;
            ok = 0;
        }
        if(nref > 1) {
            ++multref;
            ok = 0;
        }
        if(nalt > 1) {
            ++multalt;
            ok = 0;
        }
        if(naa > 1) {
            ++multaa;
            ok = 0;
        }
        if(!ok) {
            ++nbad;
            continue;
        }
#if 0
        fprintf(stderr, "  chr=%s pos=%s ref=%s alt=%s aa=%s\n",
                chr ? chr : "NULL", pos ? pos : "NULL", ref ? ref : "NULL",
                alt ? alt : "NULL", aa ? aa : "NULL");
#endif

        // Skip if ref or aa are missing.
        if(aa[0] == '.' || aa[0] == '-') {
            ++missaa;
            ok = 0;
        }
        if(ref[0] == '.' || ref[0] == '-') {
            ++missref;
            ok = 0;
        }

        if(!ok) {
            ++nbad;
            continue;
        }

        char        alleles[10];
        strcpy(alleles, ref);
        strcat(alleles, alt);
        if(strlen(alleles) != 2) {
            fprintf(stderr, "%s:%5d: Error. Number of alleles is %zu.\n",
                    __FILE__, __LINE__, strlen(alleles));
            exit(EXIT_FAILURE);
        }
        char       *aaptr = strchr(alleles, aa[0]); // ptr to ancestral allele
        if(aaptr == NULL)
            continue;
        int         aai = aaptr - alleles;  // index of ancestral allele
        assert(aai == 0 || aai == 1);

        int         x = 0, n = 0;
        gtype = strsep(&next, "\t");    // field 5

        while(gtype != NULL) {
            //# gtype is a string like "0|1" or "0/1".
            switch (gtype[0]) {
            case '.':
                break;
            case '0':
                ++n;
                break;
            case '1':
                ++x;
                ++n;
                break;
            default:
                fprintf(stderr, "%s:%d: Bad genotype: %s\n",
                        __FILE__, __LINE__, gtype);
                fprintf(stderr,
                        "  chr=%s pos=%s ref=%s alt=%s aa=%s gtype=%s\n",
                        chr ? chr : "NULL", pos ? pos : "NULL",
                        ref ? ref : "NULL", alt ? alt : "NULL",
                        aa ? aa : "NULL", gtype ? gtype : "NULL");
                exit(EXIT_FAILURE);
            }

            switch (gtype[2]) {
            case '.':
                break;
            case '0':
                ++n;
                break;
            case '1':
                ++x;
                ++n;
                break;
            default:
                fprintf(stderr, "%s:%d: Bad genotype: %s\n",
                        __FILE__, __LINE__, gtype);
                fprintf(stderr,
                        "  chr=%s pos=%s ref=%s alt=%s aa=%s gtype=%s\n",
                        chr ? chr : "NULL", pos ? pos : "NULL",
                        ref ? ref : "NULL", alt ? alt : "NULL",
                        aa ? aa : "NULL", gtype ? gtype : "NULL");
                exit(EXIT_FAILURE);
            }
            gtype = strsep(&next, "\t");    // additional fields
        }

        if(n == 0) {
            ++zerogtype;
            ++nbad;
            continue;
        } else
            ++ngood;

        if(aai == 1)
            x = n - x;
        double      p = x / ((double) n);
        printf("%4s %10s %2s %2c %20.18f\n",
               chr, pos, aa, alleles[1 - aai], p);
    }
    fprintf(stderr, "daf: %d good sites; %d rejected\n", ngood, nbad);
    if(zeroref)
        fprintf(stderr, "daf: bad sites with 0 ref alleles: %d\n", zeroref);
    if(zeroalt)
        fprintf(stderr, "daf: bad sites with 0 alt alleles: %d\n", zeroalt);
    if(zeroaa)
        fprintf(stderr, "daf: bad sites with 0 ancestral alleles: %d\n",
                zeroaa);
    if(zerogtype)
        fprintf(stderr, "daf: bad sites with 0 genotypes: %d\n", zerogtype);
    if(multref)
        fprintf(stderr, "daf: bad sites with multiple ref alleles: %d\n",
                multref);
    if(multalt)
        fprintf(stderr, "daf: bad sites with multiple alt alleles: %d\n",
                multalt);
    if(multaa)
        fprintf(stderr,
                "daf: bad sites with multiple ancestral alleles: %d\n",
                multaa);
    if(missref)
        fprintf(stderr, "daf: bad sites with missing ref alleles: %d\n",
                missref);
    if(missaa)
        fprintf(stderr, "daf: bad sites with missing ancestral alleles: %d\n",
                missaa);
    return 0;
}
