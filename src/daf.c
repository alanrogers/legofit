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

This can be generated from a vcf file that includes annotations for
ancestral alleles. If the ancestral is labelled "AA", the input for
daf can be generated, using bcftools, as follows:

  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA[\t%GT]\n' fname.vcf.gz

Output is in 5 columns, separated by whitespace:

1. chromosome
2. position of the nucleotide
3. aa, the ancestral allele
4. da, the derived allele
5. daf, derived allele frequency

The input files should include all sites at which derived alleles are
present in any of the populations under study. For example, consider
an analysis involving modern humans and Neanderthals. The modern human
data must include all sites at which Neanderthals carry derived
alleles, even if these sites do not vary among modern humans. To
accomplish this, it is best to use whole-genome data for all
populations.

The input should not contain duplicate nucleotide sites, the
chromosomes should be sorted in lexical order, and within each
chromosome, the nucleotides should be in numerical order. Otherwise,
daf will abort with an error.

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

    const int   buffsize = 16384;
    char        buff[buffsize];

    // Keep track of the number of sites at which the number
    // of reference, alternate, or ancestral alleles differs from 0
    long int         zeroref = 0, zeroalt = 0, zeroaa = 0, zerogtype = 0;
    long int         missref = 0, missaa = 0;
    long int         multref = 0, multalt = 0, multaa = 0;
    long int         indelref=0, indelalt=0, indelaa=0;
    long int         nbad = 0, ngood = 0;
    int         ok;             // is current line acceptable
    long unsigned lastnucpos = 0, nucpos;
    char        lastchr[100] = { '\0' };

    printf("#%3s %10s %2s %2s %s\n", "chr", "pos", "aa", "da", "daf");
    while(1) {
        if(NULL == fgets(buff, buffsize, stdin)) {
            break;
        }
        if(NULL == strchr(buff, '\n') && !feof(stdin)) {
            fprintf(stderr, "%s:%d: Buffer overflow. size=%zu\n",
                    __FILE__, __LINE__, sizeof(buff));
            exit(EXIT_FAILURE);
        }
        char       *chr, *pos, *reftoken, *alttoken, *aatoken;
        char       *gtype, *next = buff;

        chr = strsep(&next, "\t");  // field 0
        pos = strsep(&next, "\t");  // field 1
        reftoken = strsep(&next, "\t");  // field 2
        alttoken = strsep(&next, "\t");  // field 3
        aatoken = strsep(&next, "\t");   // field 4

        nucpos = strtoul(pos, NULL, 10);

        if(aatoken == NULL) {
            fprintf(stderr, "%s:%d: Bad input line\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        // lowercase alleles
        strlowercase(reftoken);
        strlowercase(alttoken);
        strlowercase(aatoken);

        // strip space characters
        (void) stripchr(chr, ' ');
        (void) stripchr(pos, ' ');
        (void) stripchr(reftoken, ' ');
        (void) stripchr(alttoken, ' ');
        (void) stripchr(aatoken, ' ');

        // Check sort of chromosomes
        if(*lastchr) {
            int         diff = strcmp(lastchr, chr);
            if(diff > 0) {
                // bad sort
                fprintf(stderr, "%s:%d: unsorted chromosomes\n",
                        __FILE__, __LINE__);
                fprintf(stderr, "    %s > %s\n", lastchr, chr);
                exit(1);
            } else if(diff < 0) {
                // new chromosome
                int         status =
                    snprintf(lastchr, sizeof lastchr, "%s", chr);
                if(status >= sizeof lastchr) {
                    fprintf(stderr, "%s:%d: buffer overflow\n",
                            __FILE__, __LINE__);
                    exit(1);
                }
                lastnucpos = 0;
            }
        } else {
            // initialize lastchr
            int         status = snprintf(lastchr, sizeof lastchr, "%s", chr);
            if(status >= sizeof lastchr) {
                fprintf(stderr, "%s:%d: buffer overflow\n",
                        __FILE__, __LINE__);
                exit(1);
            }
            assert(lastnucpos == 0);
        }

        // Check sort of nucleotide positions
        if(lastnucpos) {
            if(lastnucpos == nucpos) {
                fprintf(stderr, "%s:%d: Duplicate: chr=%s pos=%lu\n",
                        __FILE__, __LINE__, chr, nucpos);
                fprintf(stderr, "%s:%d: Previous : chr=%s pos=%lu\n",
                        __FILE__, __LINE__, lastchr, lastnucpos);
                exit(1);
            } else if(lastnucpos > nucpos) {
                fprintf(stderr, "%s:%d: Missorted nucleotide positions\n",
                        __FILE__, __LINE__);
                fprintf(stderr, "   Current : chr=%s pos=%lu\n", chr, nucpos);
                fprintf(stderr, "   Previous: chr=%s pos=%lu\n",
                        lastchr, lastnucpos);
                exit(1);
            }
        }
        lastnucpos = nucpos;

        int maxalleles = 10, nref, nalt, naa;
        char *ref[maxalleles];
        char *alt[maxalleles];
        char *aa[maxalleles];

        // Replace '-' with '.' in allele tokens, so there is only
        // one kind of missing value.
        strReplaceChr(reftoken, '-', '.');
        strReplaceChr(alttoken, '-', '.');
        strReplaceChr(aatoken, '-', '.');

        // REF and ALT alleles are separated by commas.
        // Ancestral alleles are separated by "|".
        nref = tokenize(maxalleles, ref, reftoken, ",");
        nalt = tokenize(maxalleles, alt, alttoken, ",");
        naa = tokenize(maxalleles, aa, aatoken, "|");

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

        // Skip if ref, alt or aa are indels
        if(strlen(ref[0]) != 1) {
            ++indelref;
            ok=0;
        }
        if(strlen(alt[0]) != 1) {
            ++indelalt;
            ok=0;
        }
        if(strlen(aa[0]) != 1) {
            ++indelaa;
            ok=0;
        }
        if(!ok) {
            ++nbad;
            continue;
        }

        // Skip if ref or aa are missing.
        if(0==strcmp(aa[0], ".")) {
            ++missaa;
            ok = 0;
        }
        if(0==strcmp(ref[0], ".")) {
            ++missref;
            ok = 0;
        }

        if(!ok) {
            ++nbad;
            continue;
        }

        // If alt is missing and aa differs from ref, then set
        // alt equal to aa.
        if(0==strcmp(alt[0], ".") && 0!=strcmp(ref[0], aa[0]))
            alt[0] = aa[0];

        char        *alleles[2];
        int          aai;
        alleles[0] = ref[0];
        alleles[1] = alt[0];

        // set index of ancestral allele
        if(0==strcmp(alleles[0], aa[0]))
            aai=0;
        else if(0==strcmp(alleles[1], aa[0]))
            aai=1;
        else {
            // skip site: there are 3 alleles
            continue;
        }
        assert(aai==0 || aai==1);

        int         x = 0, n = 0;
        gtype = strsep(&next, "\t");    // field 5

        while(gtype != NULL) {
            if(strlen(gtype != 3)) {
                fprintf(stderr,"%s:%d: Bad genotype: %s\n",
                        __FILE__, __LINE__, gtype);
                fprintf(stderr,
                        "  chr=%s pos=%s ref=%s alt=%s gtype=%s\n",
                        chr, pos, ref[0], alt[0], gtype);
                exit(EXIT_FAILURE);
            }
            // gtype is a string like "0|1" or "0/1".
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
                        chr, pos, ref[0], alt[0], aa[0], gtype);
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
                        chr, pos, ref[0], alt[0], aa[0], gtype);
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
        printf("%4s %10s %2s %2s %0.18g\n",
               chr, pos, aa[0], alleles[1 - aai], p);
    }
    fprintf(stderr, "daf: %ld good sites; %ld rejected\n", ngood, nbad);
    if(zeroref)
        fprintf(stderr, "daf: bad sites with 0 ref alleles: %ld\n", zeroref);
    if(zeroalt)
        fprintf(stderr, "daf: bad sites with 0 alt alleles: %ld\n", zeroalt);
    if(zeroaa)
        fprintf(stderr, "daf: bad sites with 0 ancestral alleles: %ld\n",
                zeroaa);
    if(zerogtype)
        fprintf(stderr, "daf: bad sites with 0 genotypes: %ld\n", zerogtype);
    if(multref)
        fprintf(stderr, "daf: bad sites with multiple ref alleles: %ld\n",
                multref);
    if(multalt)
        fprintf(stderr, "daf: bad sites with multiple alt alleles: %ld\n",
                multalt);
    if(multaa)
        fprintf(stderr,
                "daf: bad sites with multiple ancestral alleles: %ld\n",
                multaa);
    if(indelref)
        fprintf(stderr,
                "daf: bad sites with ref allele an indel: %ld\n",
                indelref);
    if(indelalt)
        fprintf(stderr,
                "daf: bad sites with alt allele an indel: %ld\n",
                indelalt);
    if(indelaa)
        fprintf(stderr,
                "daf: bad sites with ancestral allele an indel: %ld\n",
                indelaa);
    if(missref)
        fprintf(stderr, "daf: bad sites with missing ref alleles: %ld\n",
                missref);
    if(missaa)
        fprintf(stderr, "daf: bad sites with missing ancestral alleles: %ld\n",
                missaa);
    return 0;
}
