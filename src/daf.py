#!/usr/bin/python

# Calculate derived allele frequency, daf.
#
# Input file should consist of comma-separated columns:
# Col 1: chromosome
# Col 2: position
# Col 3: reference allele
# Col 4: alternate alleles
# Col 5: ancestral allele
# Cols 6..: genotypes in format "0/1" or "0|1", where 0 represents
#           a copy of the reference allele and 1 a copy of the derived
#           allele.
#
# With 1000-genomes data, this input can be generated from a vcf or
# bcf file as follows:
#    bcftools query -f '%CHROM,%POS,%REF,%ALT,%INFO/AA[,%GT]\n' fname.bcf
#
# Output is in 3 columns, separated by a space character:
# Col 1: chromosome
# Col 2: position of the nucleotide
# Col 3: aa, the ancestral allele
# Col 4: da, the derived allele
# Col 5: daf, derived allele frequency
#
# If ref, alt, or the ancestral allele consists of more than a single
# character, the site is skipped. 

import sys

print "#%3s %10s %2s %2s %20s" % ("chr", "pos", "aa", "da", "daf")
for line in sys.stdin:
    line = line.strip().lower().split(",")
    chr = line[0]
    pos = line[1]
    ref = line[2]
    alt = line[3]
    aa = line[4].replace("|","")
    if len(ref)!=1 or len(alt)!=1 or len(aa)!=1:
        continue
    alleles = ref+alt
    if len(alleles) != 2:
        print "Error number of alleles is", len(alleles)
        sys.exit(1)
    aai = alleles.find(aa)
    if aai < 0 or aai >= 2:
        # ancestral allele is neither ref nor alt.
        continue
    x = 0
    n = 0
    for gtype in line[5:]:
        # gtype is a string like "0|1" or "0/1".
        g = gtype[0]
        if g=='0' or g=='1':
            n += 1
            x += int(g)

        g = gtype[2]
        if g=='0' or g=='1':
            n += 1
            x += int(g)

    if n == 0:
        continue

    if aai == 1:
        x = n-x
    p = float(x)/float(n)

#    print line
    print "%4s %10s %2s %2s %20.18f" % (chr, pos, aa, alleles[1-aai], p)
