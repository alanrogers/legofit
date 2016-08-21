// Calculate derived allele frequency, daf.
//
// Input file should consist of comma-separated columns:
// Col 1: chromosome
// Col 2: position
// Col 3: reference allele
// Col 4: alternate alleles
// Col 5: ancestral allele
// Cols 6..: genotypes in format "0/1" or "0|1", where 0 represents
//           a copy of the reference allele and 1 a copy of the derived
//           allele.
// With 1000-genomes data, this input can be generated from a vcf or
// bcf file as follows:
// bcftools query -f '%CHROM,%POS,%REF,%ALT,%INFO/AA[,%GT]\n' fname.bcf//
//
// Output is in 3 columns, separated by a space character:
// Col 1: chromosome
// Col 2: position of the nucleotide
// Col 3: aa, the ancestral allele
// Col 4: da, the derived allele
// Col 5: daf, derived allele frequency
//
// If ref, alt, or the ancestral allele consists of more than a single
// character, the site is skipped. 
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int main(int argc, char **argv) {
    
    const int buffsize = 4096;
    char buff[buffsize];

    printf("#%3s %10s %2s %2s %20s\n", "chr", "pos", "aa", "da", "daf");
    while(1) {
        if(NULL==fgets(buff, buffsize, stdin))
            break;
        char *chr, *pos, *ref, *alt, *aa, *gtype, *next = buff;

        chr = strsep(&next, ","); // field 0
        pos = strsep(&next, ","); // field 1
        ref = strsep(&next, ","); // field 2
        alt = strsep(&next, ","); // field 3
        aa = strsep(&next, ",");  // field 4

        if(aa==NULL) {
            fprintf(stderr,"%s:%d: Bad input line\n",__FILE__,__LINE__);
            exit(EXIT_FAILURE);
        }

        // lowercase alleles
        strlowercase(ref);
        strlowercase(alt);
        strlowercase(aa);

		// strip extraneous characters
		(void) stripchr(chr, ' ');
		(void) stripchr(pos, ' ');
		int nref = stripchr(ref, ' '); // nref is number of ref alleles
		int nalt = stripchr(alt, ' '); // nalt is number of alt alleles
		(void) stripchr(aa, ' ');
        int naa = stripchr(aa, '|');   // naa is number of ancestral alleles
        if(nref != 1 || nalt != 1 || naa != 1)
            continue;

		if(aa[0] == '.' || aa[0] == '-')
			continue;

#if 0
        fprintf(stderr,"  chr=%s pos=%s ref=%s alt=%s aa=%s\n",
                chr?chr:"NULL", pos?pos:"NULL", ref?ref:"NULL",
                alt?alt:"NULL", aa?aa:"NULL");
#endif

        char alleles[10];
        strcpy(alleles, ref);
        strcat(alleles, alt);
        if(strlen(alleles) != 2) {
            fprintf(stderr, "%s:%5d: Error. Number of alleles is %zu.\n",
                    __FILE__,__LINE__,strlen(alleles));
            exit(EXIT_FAILURE);
        }
        char *aaptr = strchr(alleles, aa[0]); // ptr to ancestral allele
        if(aaptr==NULL)
            continue;
        int aai = aaptr - alleles;  // index of ancestral allele
        assert(aai==0 || aai==1);

        int x = 0, n = 0;
        gtype = strsep(&next, ",");  // field 5

        while(gtype != NULL) {
            //# gtype is a string like "0|1" or "0/1".
			switch(gtype[0]) {
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
				fprintf(stderr,"%s:%d: Bad genotype: %s\n",
						__FILE__,__LINE__, gtype);
                fprintf(stderr,"  chr=%s pos=%s ref=%s alt=%s aa=%s gtype=%s\n",
                        chr?chr:"NULL", pos?pos:"NULL", ref?ref:"NULL",
                        alt?alt:"NULL", aa?aa:"NULL", gtype?gtype:"NULL");
				exit(EXIT_FAILURE);
			}

			switch(gtype[2]) {
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
				fprintf(stderr,"%s:%d: Bad genotype: %s\n",
						__FILE__,__LINE__, gtype);
				exit(EXIT_FAILURE);
			}
			gtype = strsep(&next, ",");  // additional fields
		}

		if(n == 0)
			continue;

		if( aai == 1 )
			x = n-x;
		double p = x/((double) n);
		printf("%4s %10s %2s %2c %20.18f\n",
			   chr, pos, aa, alleles[1-aai], p);
	}
	return 0;
}
