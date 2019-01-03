/**
 * @file boot.c
 * @author Alan R. Rogers
 * @brief Functions for a moving blocks bootstrap.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "boot.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

/// Contains the all data involved in a moving blocks bootstrap
struct Boot {
    int nchr;              ///< number of chromosomes
    long nsnp;             ///< number of snps (across all chromosomes)
    long blocksize;        ///< number of SNPs per block
    long nrep;             ///< number of bootstrap replicates 
    long nblock;           ///< number of blocks
    int npat;              ///< number of site patterns
    double **count;        ///< count[i][j]: j'th site pattern in i'th rep
    long **start;          ///< start[i][j] = start of j'th block in i'th rep

    /// cum[i] is number of SNPs preceding chromosome i
    long *cum;
};

/// Contains the data for a bootstrap confidence interval.
struct BootConf {
    long nrep;                  ///< repetitions
    long blocksize;             ///< nucleotide positions per block
    double confidence;          ///< size of confidence region
    double *low, *high;         ///< confidence bounds
};

long   LInt_div_round(long num, long denom);
double interpolate(double p, double *v, long len);
long   adjustBlockLength(long lengthWanted, int nsnp);
long   Boot_multiplicity(const Boot * self, long snpndx, long rep);

/// Divide num by denom and round the result to the nearest integer.
long LInt_div_round(long num, long denom) {
    assert(denom != 0L);
    ldiv_t quotrem = ldiv(num, denom);
    if(2L * quotrem.rem > denom)
        return 1L + quotrem.quot;
    return quotrem.quot;
}

/// Return a blocksize that is as close as possible to lengthWanted
/// while still making length*nblock close to nsnp.
long adjustBlockLength(long lengthWanted, int nsnp) {
    long nblock = LInt_div_round(nsnp, lengthWanted);
    return LInt_div_round(nsnp, nblock);
}

/// Constructor for class Boot.
Boot *Boot_new(int nchr, long nsnpvec[nchr], long nrep, int npat,
                     long blocksize, gsl_rng * rng) {
    long i, j;
    assert(blocksize > 0);
    if(nrep == 0)
        return NULL;

    Boot *self = malloc(sizeof(Boot));
    CHECKMEM(self);
    self->nchr = nchr;

    self->cum = malloc(nchr * sizeof(self->cum[0]));
    CHECKMEM(self->cum);

    // self->nsnp is the total number of SNPs, summed across chromosomes
    self->nsnp = 0;
    for(i=0; i < nchr; ++i)
        self->nsnp += nsnpvec[i];

    // self->cum[i] is number of SNPs preceding chr i.
    self->cum[0] = 0;
    for(i=1; i < nchr; ++i)
        self->cum[i] = self->cum[i-1] + nsnpvec[i-1];

    self->nrep = nrep;
    self->npat = npat;
    self->blocksize = blocksize = adjustBlockLength(blocksize, self->nsnp);

    if(self->blocksize > self->nsnp) {
        fprintf(stderr,
                "%s:%s:%d: blocksize must be < nsnp.\n"
                "    However, blocksize=%ld and nsnp=%ld\n",
                __FILE__, __func__, __LINE__,
                self->blocksize, self->nsnp);
        fprintf(stderr, " Use --blocksize argument"
                " to reduce blocksize.\n");
        exit(EXIT_FAILURE);
    }

    self->nblock = LInt_div_round(self->nsnp, self->blocksize);

    // Block start positions are uniform on [0, nsnp-blocksize+1).
    unsigned long endpos = self->nsnp - self->blocksize + 1;

    // Allocate arrays
    self->start = calloc((unsigned long) self->nrep, sizeof(self->start[0]));
    CHECKMEM(self->start);

    self->count = calloc((unsigned long) self->nrep, sizeof(self->count[0]));
    CHECKMEM(self->count);

    for(i = 0; i < self->nrep; ++i) {
        self->start[i] = calloc((unsigned long) self->nblock,
                                sizeof(self->start[0][0]));
        CHECKMEM(self->start[i]);

        self->count[i] = calloc((unsigned long) self->npat,
                                sizeof(self->count[i][0]));
        CHECKMEM(self->count[i]);
    }

    // For each replicate, define sorted array of block starting positions.
    for(i = 0; i < self->nrep; ++i) {
        for(j = 0; j < self->nblock; ++j)
            self->start[i][j] = gsl_rng_uniform_int(rng, endpos);

        qsort(self->start[i], (size_t) self->nblock,
              sizeof(self->start[0][0]), compareLongs);
    }

#ifndef NDEBUG
    Boot_sanityCheck(self, __FILE__, __LINE__);
#endif
    return self;
}

#ifndef NDEBUG
void Boot_sanityCheck(const Boot * self, const char *file, int line) {
    long i, j;
    REQUIRE(self->nchr > 0, file, line);
    REQUIRE(self->nsnp >= self->cum[self->cum[self->nchr - 1]],
            file, line);
    REQUIRE(self->cum[0] == 0, file, line);
    for(i=1; i < self->nchr; ++i)
        REQUIRE(self->cum[i] > self->cum[i-1], file, line);
    REQUIRE(self->blocksize > 0, file, line);
    REQUIRE(self->blocksize < 100000, file, line);
    REQUIRE(self != NULL, file, line);
    REQUIRE(self->nrep > 0, file, line);
    REQUIRE(self->nsnp > 0, file, line);
    REQUIRE(self->nblock > 0, file, line);
    REQUIRE(self->npat > 0, file, line);

    REQUIRE(self->count != NULL, file, line);
    REQUIRE(self->start != NULL, file, line);

    unsigned long endpos = self->nsnp - self->blocksize + 1;
    long prev;

    for(i = 0; i < self->nrep; ++i) {
        REQUIRE(self->count[i] != NULL, file, line);
        REQUIRE(self->start[i] != NULL, file, line);
        for(j = 0; j < self->nblock; ++j) {
            REQUIRE(self->start[i][j] >= 0, file, line);
            REQUIRE(self->start[i][j] < endpos, file, line);
            if(j > 0)
                REQUIRE(self->start[i][j] >= prev, file, line);
            prev = self->start[i][j];
        }
        for(j = 0; j < self->npat; ++j)
            REQUIRE(self->count[i][j] >= 0.0, file, line);
    }
}
#endif

/// How many copies of snp with index snpndx are present in a given
/// repetition (rep)? 
long Boot_multiplicity(const Boot * self, long snpndx, long rep) {
    long lndx, hndx, lowtarget;

    assert(snpndx < self->nsnp);

    // lndx is index of first block containing snp
    lowtarget = snpndx - self->blocksize + 1;
    lndx = long_first_geq(lowtarget, self->start[rep], self->nblock);
    if(lndx == self->nblock || self->start[rep][lndx] > snpndx)
        return 0;

    assert(snpndx >= self->start[rep][lndx]);
    assert(snpndx - self->start[rep][lndx] < self->blocksize);

    // hndx is index of first block not containing snp
    // First line below searches the sub-array beginning with
    // entry lndx. This returns an index into the sub-array.
    // The second line adds lndx to generate an index into the full
    // array.  
    hndx = long_first_geq(snpndx + 1, self->start[rep] + lndx,
                          self->nblock - lndx);
    hndx += lndx;

    assert(hndx == 0
           || self->start[rep][hndx - 1] - snpndx < self->blocksize);

    return hndx - lndx;
}

/**
 * Add one site pattern contribution to a Boot structure.
 * @param [inout] self The Boot structure to modify.
 * @param [in] chr The index of the chromosome to modify.
 * @param [in] snpndx The index of the current snp.
 * @param [in] pat The index of the current site pattern.
 * @param [in] z the contribution of the snp to the site pattern.
 */
void Boot_add(Boot * self, int chr, long snpndx, int pat, double z) {
    assert(pat < self->npat);
    assert(snpndx < self->nsnp);
    assert(z >= 0.0);
    snpndx  += self->cum[chr];
    for(register int rep = 0; rep < self->nrep; ++rep) {

        // w is the number times the current snp is represented
        // in the current bootstrap replicate.
        register long w = Boot_multiplicity(self, snpndx, rep);

        self->count[rep][pat] += w * z;
    }
}

/// Destructor
void Boot_free(Boot * self) {
#ifndef NDEBUG
    Boot_sanityCheck(self, __FILE__, __LINE__);
#endif

    for(int i = 0; i < self->nrep; ++i) {
        free(self->start[i]);
        free(self->count[i]);
    }
    free(self->start);
    free(self->count);
    free(self->cum);
    free(self);
}

/// Add to an array the site pattern counts from a bootstrap
/// replicate.
/// @param [in] self Points to a Boot object.
/// @param [in] the index of the bootstrap replicate
/// @param [in] npat the number of site patterns
/// @param [out] count An array of doubles. The function will add to
/// count[i] the contribution of site pattern i in bootstrap replicate
/// rep. 
void Boot_aggregate(Boot * self, int rep, int npat, double count[npat]) {
    assert(self);
    assert(npat == self->npat);
    int i;

#ifndef NDEBUG
    for(i = 0; i < npat; ++i) {
        if(count[i] != 0.0) {
            fprintf(stderr, "%s:%s:%d: count not initialized.\n",
                    __FILE__, __func__, __LINE__);
            exit(EXIT_FAILURE);
        }
    }
#endif

    for(i = 0; i < self->npat; ++i)
        count[i] += self->count[rep][i];
}

/// Interpolate in order to approximate the value v[p*(len-1)].
/// Return NaN if len==0.
double interpolate(double p, double *v, long len) {
    if(len == 0)
        return strtod("NAN", 0);
    long i, j;
    double w;
    double goal = p * (len - 1);

    i = floor(goal);
    j = ceil(goal);

    assert(i >= 0);
    assert(j < len);

    if(i == j)                  /* no interpolation needed */
        return v[i];
    w = goal - i;
    return (1.0 - w) * v[i] + w * v[j];
}

/**
 * Calculate confidence bounds from a vector of values representing
 * samples drawn from the sampling distribution of some estimator.
 *
 * To calculate the lower bound (*lowBnd), the function calculates the
 * total probability mass in the tails (1 - confidence) and divides
 * this into two equal parts to find p, the probability mass in each
 * tail. It then estimates a value L such that a fraction p of the data
 * values are less than or equal to L. To find this value, the
 * function uses linear interpolation between the sorted list of data
 * values.
 *
 * The upper bound (*highBnd) is calculated in an analogous fashion.
 *
 * @param[out] lowBnd,highBnd Calculated results will be written into
 * these memory locations.
 * @param[in] confidence Fraction of sampling distribution that lies
 * inside the confidence bounds.
 * @param[in] len The number of values inf v.
 * @param[in] v The vector of values.
 * @sideeffect Sorts the vector v.
 */
void confidenceBounds(double *lowBnd, double *highBnd, double confidence,
                      long len, double v[len]) {
    double tailProb = (1.0 - confidence) / 2.0;

    qsort(v, (size_t) len, sizeof(v[0]), compareDoubles);
    *lowBnd = interpolate(tailProb, v, len);
    *highBnd = interpolate(1.0 - tailProb, v, len);
}

/// Print a Boot object
void Boot_print(const Boot * self, FILE * ofp) {
    long rep, j;

    fprintf(ofp,
            "Boot_print: nsnp=%ld nrep=%ld blocksize=%ld nblock=%ld\n",
            self->nsnp, self->nrep, self->blocksize, self->nblock);

    fprintf(ofp, "Block starts:\n");
    for(rep = 0; rep < self->nrep; ++rep) {
        fprintf(ofp, "  rep %ld:", rep);
        for(j = 0; j < self->nblock; ++j)
            fprintf(ofp, " %ld", self->start[rep][j]);
        putc('\n', ofp);
    }

    fprintf(ofp, "Site pattern counts:\n");
    for(rep = 0; rep < self->nrep; ++rep) {
        fprintf(ofp, "  rep %ld:", rep);
        for(j = 0; j < self->npat; ++j)
            fprintf(ofp, " %lf", self->count[rep][j]);
        putc('\n', ofp);
    }
}

#ifndef NDEBUG

unsigned Boot_multiplicity_slow(Boot * self, long snp, long rep);

/** For debugging Boot_multiplicity */
unsigned Boot_multiplicity_slow(Boot * self, long snp, long rep) {
    unsigned i, n = 0;

    for(i = 0; i < self->nblock; ++i) {
        long distance = snp - self->start[rep][i];

        if(distance < 0)
            break;
        if(distance < self->blocksize)
            ++n;
    }
    return n;
}
#endif
