/**
@file raf.c
@page raf
@brief Calculate reference allele frequency, raf.

Input file should consist of tab-separated columns:

1. chromosome
2. position
3. reference allele
4. alternate alleles
5. genotype in format "0/1" or "0|1", where 0 represents
a copy of the reference allele and 1 a copy of the derived
allele.
6. etc for as many columns as there are genotypes.

This can be generated from a vcf file as follows:

  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' fname.vcf.gz

Output is in 5 columns, separated by tabs:

1. chromosome
2. position of the nucleotide
3. ref, the reference allele
4. alt, the alternate allele
5. raf, reference allele frequency

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
raf will abort with an error.

Sites are rejected unless they have a single ref. Missing values are
allowed for the alt allele. At the end of the job a summary of
rejected sites is written to stderr.

@copyright Copyright (c) 2017, Alan R. Rogers
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
        fprintf(stderr, "Usage: raf\n");
        fprintf(stderr, "       Reads standard input;"
                " writes to standard output.\n");
        exit(EXIT_FAILURE);
    }

    const int   buffsize = 16384;
    char        buff[buffsize];

    // Keep track of the number of sites at which the number
    // of reference or alternate alleles differs from 0
    int         zeroref = 0, zeroalt = 0, zerogtype = 0;
    int         missref = 0;
    int         multref = 0, multalt = 0;
    int         nbad = 0, ngood = 0;
    int         ok;             // is current line acceptable
    long unsigned lastnucpos = 0, nucpos;
    char        lastchr[100] = { '\0' };

    printf("#%s\t%s\t%s\t%s\t%s\n", "chr", "pos", "ref", "alt", "raf");
    while(1) {
        if(NULL == fgets(buff, buffsize, stdin)) {
            break;
        }
        if(NULL == strchr(buff, '\n') && !feof(stdin)) {
            fprintf(stderr, "%s:%d: Buffer overflow. size=%zu\n",
                    __FILE__, __LINE__, sizeof(buff));
            exit(EXIT_FAILURE);
        }
        char       *chr, *pos, *reftoken, *alttoken;
        char       *gtype, *next = buff;

        chr = strsep(&next, "\t");  // field 0
        pos = strsep(&next, "\t");  // field 1
        reftoken = strsep(&next, "\t");  // field 2
        alttoken = strsep(&next, "\t");  // field 3

        nucpos = strtoul(pos, NULL, 10);

        // lowercase alleles
        strlowercase(reftoken);
        strlowercase(alttoken);

        // strip space characters
        (void) stripchr(chr, ' ');
        (void) stripchr(pos, ' ');
        (void) stripchr(reftoken, ' ');
        (void) stripchr(alttoken, ' ');

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

        int maxalleles = 10, nref, nalt;
        char *ref[maxalleles];
        char *alt[maxalleles];

        // Replace '-' with '.' in allele tokens, so there is only
        // one kind of missing value.
        strReplaceChr(reftoken, '-', '.');
        strReplaceChr(alttoken, '-', '.');

        // REF and ALT alleles are separated by commas.
        nref = tokenize(maxalleles, ref, reftoken, ",");
        nalt = tokenize(maxalleles, alt, alttoken, ",");

        // Skip sites at which the number of reference or alternate
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
        if(nref > 1) {
            ++multref;
            ok = 0;
        }
        if(nalt > 1) {
            ++multalt;
            ok = 0;
        }
        if(!ok) {
            ++nbad;
            continue;
        }
        // Skip if ref is missing.
        if(0==strcmp(ref[0], ".")) {
            ++missref;
            ok = 0;
        }

        if(!ok) {
            ++nbad;
            continue;
        }

        int         x = 0, n = 0;
        gtype = strsep(&next, "\t");    // field 4

        while(gtype != NULL) {
            //# gtype is a string like "0|1" or "0/1".
            switch (gtype[0]) {
            case '.':
                break;
            case '0':
                ++x;
                ++n;
                break;
            case '1':
                ++n;
                break;
            default:
                fprintf(stderr, "%s:%d: Bad genotype: %s\n",
                        __FILE__, __LINE__, gtype);
                fprintf(stderr,
                        "  chr=%s pos=%s ref=%s alt=%s gtype=%s\n",
                        chr, pos, ref[0], alt[0], gtype);
                exit(EXIT_FAILURE);
            }

            switch (gtype[2]) {
            case '.':
                break;
            case '0':
                ++x;
                ++n;
                break;
            case '1':
                ++n;
                break;
            default:
                fprintf(stderr, "%s:%d: Bad genotype: %s\n",
                        __FILE__, __LINE__, gtype);
                fprintf(stderr,
                        "  chr=%s pos=%s ref=%s alt=%s gtype=%s\n",
                        chr, pos, ref[0], alt[0], gtype);
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

        double      p = x / ((double) n);
        printf("%s\t%s\t%s\t%s\t%0.18g\n",
               chr, pos, ref[0], alt[0], p);
    }
    fprintf(stderr, "raf: %d good sites; %d rejected\n", ngood, nbad);
    if(zeroref)
        fprintf(stderr, "raf: bad sites with 0 ref alleles: %d\n", zeroref);
    if(zeroalt)
        fprintf(stderr, "raf: bad sites with 0 alt alleles: %d\n", zeroalt);
    if(zerogtype)
        fprintf(stderr, "raf: bad sites with 0 genotypes: %d\n", zerogtype);
    if(multref)
        fprintf(stderr, "raf: bad sites with multiple ref alleles: %d\n",
                multref);
    if(multalt)
        fprintf(stderr, "raf: bad sites with multiple alt alleles: %d\n",
                multalt);
    if(missref)
        fprintf(stderr, "raf: bad sites with missing ref alleles: %d\n",
                missref);
    return 0;
}
