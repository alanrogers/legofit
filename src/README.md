/** 
\mainpage 

# Introduction 

Lego is a computer package that uses counts of nucleotide site
patterns to estimate the history of population size, subdivision, and
gene flow.

## Site patterns, counts, and their expectations

Consider a sample consisting of one haploid genome drawn from each of
3 populations, *X*, *Y*, and *Z*. Suppose that, at a given nucleotide
site, the derived allele is present in the genomes from *X* and *Y*
but not that from *Z*. If so, then this nucleotide position will be
said to exhibit the "*xy* site pattern." We ignore cases in which the
derived allele is present in none of the samples, in only one of them,
or in all of them. In other words, we consider only polymorphic,
non-singleton site patterns. For the special case of the 3-population
sample just described, there are only 3 such site patterns: *xy*,
*xz*, and *yz*.

In the general case, with samples from *K* populations, the number of
site patterns is \f$2^K - K - 2\f$. For example, there are 10 site
patterns in a sample involving \f$K=4\f$ populations. The table below
shows data from a sample involving 4 populations, *X*, *Y*, *N*, and
*D*. 

    SitePat          E[count]
         xy    340952.4592501
         xn     46874.1307236
         xd     46034.4670204
         yn     55137.4236715
         yd     43535.5248078
         nd    231953.3372578
        xyn     91646.1277991
        xyd     88476.9619569
        xnd     96676.3877423
        ynd    100311.4411513

The `E[count]` column shows numbers that can be thought of loosely as
counts of site patterns in a genome-wide sample. The last line tells
us that the *ynd* site pattern occurs at over 100,000 nucleotide
sites.

These number cannot really be counts, because they aren't
integers. This reflects the fact that our sample includes more than
one haploid genome from each population, and a given SNP may
contribute to several site patterns. The contribution to a given site
pattern is the probability that a sub-sample, consisting of one
haploid genome drawn at random from the larger sample of each
population, would exhibit this site pattern. For example, suppose we
have samples from three populations, *X*, *Y*, and *N*, and let
\f$p_{iX}\f$, \f$p_{iY}\f$, and \f$p_{iN}\f$ represent the frequencies
of the derived allele at the *i*th SNP in these three
samples. Then site pattern *xy* occurs at SNP *i* with
probability \f$z_i = p_{iX}p_{iY}(1-p_{iN})\f$ (Patterson et al 2010,
Science, 328(5979):S129).  Aggregating over SNPs, \f$I_{xy} = \sum_i
z_i\f$ summarizes the information in the data about this site
pattern. These are the numbers that appear in the 2nd column of the
table above.

# Installation and testing

The package is available at
[github](github.com/alanrogers/lego). Before compiling, you must
install two libraries: `pthreads` and
[`gsl`](http://www.gnu.org/software/gsl). You will need not only the
libraries themselves but also several header files, such as
`pthread.h`. I didn't need to install `pthreads`, because it came
bundled with the Gnu C compiler. But the gsl was an extra. Under
ubuntu Linux, you can install it like this:

    sudo apt-get install libgsl0-dev

On the mac, using homebrew, the command is

    brew install gsl

By default, the executable files will be copied into a directory named
`bin` in your home directory. If you want them to go somewhere else,
edit the first non-comment line of src/Makefile.

Then 

1. Cd into the src directory.
2. Type "make".
3. Type "make install".

This will try to place the executables into directory "bin" in the
user's home directory. Make sure this directory appears in your
PATH, so that the shell can find it.

This installation will work under unix-like operating systems, such as
linux and Apple's osx. I haven't tried to port this software to
Windows. 

The directory `test` contains a unit test for many of the .c files in
directory `src`. Within this directory, type

1. make xboot
2. ./xboot

to test the source file `boot.c`.  To run all unit tests, type
"make". This will take awhile, as some of the unit tests are slow.

# Genetic input data

Before doing data analysis with `lego`, you must generate data files
in "daf" format. Such files end with ".daf", which stands for "derived
allele frequency. Here are the first few lines of one such file:

    #chr        pos aa da                  daf
       1     752566  g  a 0.835294117647058854
       1     754192  a  g 0.858823529411764652
       1     755225  t  g 0.000000000000000000
       1     755228  t  g 0.000000000000000000
       1     765437  g  a 0.000000000000000000

The first line (beginning with "#") is an optional comment, which is
used here to label the columns. The columns are as follows:

1. Character strings that label chromosomes or scaffolds.
2. Position of the SNP on the chromosome or scaffold, measured in base
   pairs. Daf format doesn't care whether nucleotide positions are
   numbered beginning with 0 or with 1, provided that they are consistent
   across files in a given analysis.
4. Ancestral allele, a single letter.
5. Derived allele, also a single letter. Loci with 3 or more alleles
   should be excluded.
6. Frequency of the derived allele within the sample.

The lines should be sorted lexically by chromosome. Within
chromosomes, they should be sorted in ascending numerical order of
column 2.

*/
