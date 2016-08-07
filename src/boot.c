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

/// Contains the all data involved in a moving blocks bootstrap of
/// a single chromosome.
struct BootChr {
    long        blockLength;    // number of SNPs per block
    long        nrep;           // number of bootstrap replicates 
    long        nsnp;           // number of snps
    long        nblock;         // number of blocks
    int         npat;           // number of site patterns
    double    **count;          // count[i][j]: j'th site pattern in i'th rep
    long      **start;          // start[i][j] = start of j'th block in i'th rep
};

struct Boot {
    int nchr;                   // number of chromosomes
    BootChr **bc;               // bc[i]: bootstrap for i'th chromosome
};

/** Contains the data for a bootstrap confidence interval. */
struct BootConf {
    long        nrep;           // repetitions
    long        blockLength;    // nucleotide positions per block
    double      confidence;     // size of confidence region
    double     *low, *high;     // confidence bounds
};

long LInt_div_round(long num, long denom);
static void BootChr_allocArrays(BootChr * self);

// Divide num by denom and round the result to the nearest integer.
long LInt_div_round(long num, long denom) {
        assert(denom != 0L);
        ldiv_t quotrem = ldiv(num, denom);
        if(2L * quotrem.rem > denom)
                return 1L + quotrem.quot;
        return quotrem.quot;
}

// Return a blockLength that is as close as possible to lengthWanted
// while still making length*nblock close to nsnp.
long adjustBlockLength(long lengthWanted, int nsnp) {
	long nblock = LInt_div_round(nsnp, lengthWanted);
	return LInt_div_round(nsnp, nblock);
}

/// Constructor for class BootChr.
BootChr       *BootChr_new(long nsnp, long nrep, int npat, long blockLength,
                           gsl_rng * rng) {
    long i, j;
    assert(blockLength > 0);
    if(nrep == 0)
        return NULL;

    if(blockLength > nsnp) {
        fprintf(stderr,
                "ERR@%s:%d: in BootChr_new, nsnp must be >blockLength.\n"
                " Instead, nsnp=%ld, blockLength=%ld.\n",
                __FILE__, __LINE__, nsnp, blockLength);
        fprintf(stderr, " Use --blocksize argument"
                " to reduce blockLength.\n");
        exit(1);
    }

    // Block positions are uniform on [0, nsnp-blockLength+1).
    unsigned long endpos;
    endpos = nsnp - blockLength + 1;

    BootChr       *self = malloc(sizeof(BootChr));
    CHECKMEM(self);

    self->nsnp = nsnp;
    self->nrep = nrep;
	self->blockLength = adjustBlockLength(blockLength, nsnp);
	fprintf(stderr,"%s:%d: blockLength %ld -> %ld\n",
			__FILE__,__LINE__, blockLength, self->blockLength);
    self->npat = npat;
    self->nblock = LInt_div_round(nsnp, blockLength);

    BootChr_allocArrays(self);

    for(i = 0; i < self->nrep; ++i) {
        for(j = 0; j < self->nblock; ++j)
            self->start[i][j] = gsl_rng_uniform_int(rng, endpos);

        qsort(self->start[i], (size_t) self->nblock,
              sizeof(self->start[0][0]), compareLongs);
    }

#ifndef NDEBUG
    BootChr_sanityCheck(self, __FILE__, __LINE__);
#endif
    return self;
}

/// Allocate BootChr's arrays.
static void BootChr_allocArrays(BootChr * self) {

    long        i;

    self->start = calloc((unsigned long) self->nrep, sizeof(self->start[0]));
    CHECKMEM(self->start);

    self->count = calloc((unsigned long) self->nrep, sizeof(self->count[0]));
    CHECKMEM(self->count);

    for(i = 0; i < self->nrep; ++i) {
        self->start[i] = calloc((unsigned long) self->nblock,
                                sizeof(self->start[0][0]));
        CHECKMEM(self->start[i]);

        self->count[i] = calloc((unsigned long) self->npat,
                                sizeof(self->count[0]));
        CHECKMEM(self->count[i]);
    }
}

#ifndef NDEBUG
void BootChr_sanityCheck(const BootChr * self, const char *file, int line) {
    long        i, j;
    REQUIRE(self->blockLength > 0, file, line);
    REQUIRE(self->blockLength < 100000, file, line);
    REQUIRE(self != NULL, file, line);
    REQUIRE(self->nrep > 0, file, line);
    REQUIRE(self->nsnp > 0, file, line);
    REQUIRE(self->nblock > 0, file, line);
    REQUIRE(self->npat > 0, file, line);

    REQUIRE(self->count != NULL, file, line);
    REQUIRE(self->start != NULL, file, line);

    unsigned long endpos = self->nsnp - self->blockLength + 1;
    long        prev;

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

/// How many copies of snp are present in a given repetition (rep)?
long BootChr_multiplicity(const BootChr * self, long snpndx, long rep) {
    long        lndx, hndx, lowtarget;

    assert(snpndx < self->nsnp);

    // lndx is index of first block containing snp
    lowtarget = snpndx - self->blockLength + 1;
    lndx = long_first_geq(lowtarget, self->start[rep], self->nblock);
    if(lndx == self->nblock || self->start[rep][lndx] > snpndx)
        return 0;

    assert(snpndx >= self->start[rep][lndx]);
    assert(snpndx - self->start[rep][lndx] < self->blockLength);

    // hndx is index of first block not containing snp
    hndx = long_first_geq(snpndx + 1, self->start[rep] + lndx,
                          self->nblock - lndx);
    hndx += lndx;

    assert(hndx == 0
             || self->start[rep][hndx - 1] - snpndx < self->blockLength);

    return hndx - lndx;
}

/*
 * Add one site pattern contribution to a BootChr structure. On entry,
 * self points to the BootChr structure, snpndx is the index of the
 * current snp, pat is that of the current site pattern, and z is the
 * contribution of the snp to the site pattern.
 */
void BootChr_add(BootChr * self, long snpndx, int pat, double z) {
    assert(pat < self->npat);
    assert(snpndx < self->nsnp);
	if(!(z >= 0))
		fprintf(stderr,"%s:%s:%d: z=%lf\n", __FILE__,__func__,__LINE__,z);
    assert(z >= 0.0);
    for(register int rep = 0; rep < self->nrep; ++rep) {

        // w is the number times the current snp is represented
        // in the current bootstrap replicate.
        register long w = BootChr_multiplicity(self, snpndx, rep);

        self->count[rep][pat] += w*z;
    }
}

/// Return number of bootstrap repetitions
long BootChr_nrep(const BootChr * self) {
    assert(self);
    return self->nrep;
}

/// Return number of bootstrap repetitions
long BootChr_npat(const BootChr * self) {
    assert(self);
    return self->npat;
}

/// Return number of blocks
long BootChr_nblock(const BootChr * self) {
    assert(self);
    return self->nblock;
}

/// Return number of SNPs
long BootChr_nsnp(const BootChr * self) {
    assert(self);
    return self->nsnp;
}

// Destructor
void BootChr_free(BootChr * self) {
#ifndef NDEBUG
    BootChr_sanityCheck(self, __FILE__, __LINE__);
#endif

    for(int i = 0; i < self->nrep; ++i) {
        free(self->start[i]);
        free(self->count[i]);
    }
    free(self->start);
    free(self->count);
    free(self);
}

/// Add to array "count" the values for bootstrap repetition "rep".
void BootChr_aggregate(BootChr * self, int rep, int npat, double count[npat]) {
    assert(self);
    assert(npat == self->npat);
    int j;
    for(j=0; j < self->npat; ++j)
        count[j] += self->count[rep][j];
}

Boot * Boot_new(int nchr, long nsnp[nchr], long nrep, int npat,
                long blockLength, gsl_rng *rng) {
    Boot *self = malloc(sizeof(Boot));
    CHECKMEM(self);
    self->nchr = nchr;
    self->bc = malloc(nchr * sizeof(BootChr));
    CHECKMEM(self->bc);

    for(int i=0; i < nchr; ++i) {
        self->bc[i] = BootChr_new(nsnp[i], nrep, npat, blockLength, rng);
        CHECKMEM(self->bc[i]);
    }
    return self;
}

void Boot_free(Boot *self) {
    for(int i=0; i < self->nchr; ++i)
        free(self->bc[i]);
    free(self->bc);
    free(self);
}

void Boot_add(Boot *self, int chr, long snpndx, int pat, double z) {
    BootChr_add(self->bc[chr], snpndx, pat, z);
}

void Boot_aggregate(Boot * self, int rep, int npat,
                    double count[npat]) {
    for(int i=0; i < self->nchr; ++i)
        BootChr_aggregate(self->bc[i], rep, npat, count);
}

void Boot_sanityCheck(const Boot * self, const char *file, int line) {
    for(int i=0; i < self->nchr; ++i)
        BootChr_sanityCheck(self->bc[i], file, line);
}

/// Interpolate in order to approximate the value v[p*(len-1)].
/// Return NaN if len==0.
double interpolate(double p, double *v, long len) {
    if(len == 0)
        return strtod("NAN", 0);
    long        i, j;
    double      w;
    double      goal = p * (len - 1);

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
 * @param[in] v The vector of values.
 * @param[in] len The number of values inf v.
 * @sideeffect The function sorts the vector v.
 */
void confidenceBounds(double *lowBnd, double *highBnd, double confidence,
                      double *v, long len) {
    double      tailProb = (1.0 - confidence) / 2.0;

    qsort(v, (size_t) len, sizeof(v[0]), compareDoubles);
    *lowBnd = interpolate(tailProb, v, len);
    *highBnd = interpolate(1.0 - tailProb, v, len);
}

/// Print a BootChr object
void BootChr_print(const BootChr * self, FILE * ofp) {
    long        rep, j;

    fprintf(ofp,
            "BootChr_print: nsnp=%ld nrep=%ld blockLength=%ld nblock=%ld\n",
            self->nsnp, self->nrep, self->blockLength, self->nblock);

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
/* For debugging BootChr_multiplicity */
unsigned BootChr_multiplicity_slow(BootChr * self, long snp, long rep) {
    unsigned    i, n = 0;

    for(i = 0; i < self->nblock; ++i) {
        long        distance = snp - self->start[rep][i];

        if(distance < 0)
            break;
        if(distance < self->blockLength)
            ++n;
    }
    return n;
}
#endif

#if 0
/*
 * Allocate and initialize a duplicate of a BootChr structure. Deallocate
 * with BootChr_free.
 */
BootChr       *BootChr_dup(const BootChr * old) {
    assert(old != NULL);
    assert(old->nblock > 0);
    BootChr       *new = memdup(old, sizeof(*old));

    assert(new->nblock > 0);

    new->start = malloc(new->nrep * sizeof(new->start[0]));
    checkmem(new->start, __FILE__, __LINE__);

    new->tab = malloc(new->nrep * sizeof(new->tab[0]));
    checkmem(new->tab, __FILE__, __LINE__);

    new->spectab = malloc(new->nrep * sizeof(new->spectab[0]));
    checkmem(new->spectab, __FILE__, __LINE__);

    assert(sizeof(new->start[0][0]) > 0);

    for(int i = 0; i < new->nrep; ++i) {
        new->tab[i] = Tabulation_dup(old->tab[i]);
        new->spectab[i] = Spectab_dup(old->spectab[i]);
        new->start[i] = memdup(old->start[i],
                               new->nblock * sizeof(new->start[0][0]));
    }

#ifndef NDEBUG
    BootChr_sanityCheck(new, __FILE__, __LINE__);
#endif
    return new;
}

/* Return 1 if the two BootChr structs are idential; zero otherwise */
int BootChr_equals(const BootChr * x, const BootChr * y) {
    if(x == NULL || y == NULL)
        return false;
    if(x->nsnp != y->nsnp)
        return false;
    if(x->nrep != y->nrep)
        return false;
    if(x->blockLength != y->blockLength)
        return false;
    if(x->nBins != y->nBins)
        return false;

    long        i, j;

    for(i = 0; i < x->nrep; ++i) {
        if(!Tabulation_equals(x->tab[i], y->tab[i]))
            return false;

        if(!Spectab_equals(x->spectab[i], y->spectab[i]))
            return false;

        for(j = 0; j < x->nblock; ++j)
            if(x->start[i][j] != y->start[i][j])
                return false;
    }
    return true;
}

/*
 * Add the contents of BootChr structure y to BootChr structure x.
 * On return x will summarize all the LD values that were originally
 * contained in the two original BootChr structures.
 */
void BootChr_plus_equals(BootChr * x, const BootChr * y) {
    int         rep;

#ifndef NDEBUG
    BootChr_sanityCheck(x, __FILE__, __LINE__);
    BootChr_sanityCheck(y, __FILE__, __LINE__);
#endif
    if(x == NULL)
        die("BootChr_plus_equals: destination is NULL", __FILE__, __LINE__);
    if(y == NULL)
        return;

    if(x->nsnp != y->nsnp
       || x->nrep != y->nrep
       || x->blockLength != y->blockLength
       || x->nblock != y->nblock || x->nBins != y->nBins)
        die("BootChr_plus_equals: unconformable arguments", __FILE__, __LINE__);

    for(rep = 0; rep < x->nrep; ++rep) {
        long        block;

        for(block = 0; block < x->nblock; ++block) {
            if(x->start[rep][block] != y->start[rep][block])
                die("BootChr_plus_equals: incompatible block starts",
                    __FILE__, __LINE__);
        }
        Tabulation_plus_equals(x->tab[rep], y->tab[rep]);
        Spectab_plus_equals(x->spectab[rep], y->spectab[rep]);
    }
#ifndef NDEBUG
    BootChr_sanityCheck(x, __FILE__, __LINE__);
#endif
    return;
}

/*
 * Return the raw counts for a single rep and bin. Counts for
 * numerator, denominator, and sep_cm are returned in the
 * corresponding arguments. The function itself returns nobs.
 */
long unsigned BootChr_rawCounts(const BootChr * self, int rep, int bin,
                             double *numerator, double *denominator,
                             double *sumRsq, double *sep_cm) {

    assert(rep < self->nrep);
    assert(bin < self->nBins);

    return Tabulation_rawCounts(self->tab[rep], bin, numerator,
                                denominator, sumRsq, sep_cm);
}

/**
 * Remove replicates in which some bins have zero observations.
 * These would generate NaN values in the calculation of sigdsq,
 * which are not handled by the minimizer. Return revised number of
 * bootstrap replicates.
 */
long BootChr_purge(BootChr * self) {
    assert(self);

    long        rep, nGoodReps = self->nrep;
    int         bin;

    rep = 0;
    while(rep < nGoodReps) {
        int         clean = 1;

        for(bin = 0; bin < self->nBins; ++bin) {
            if(Tabulation_nObs(self->tab[rep], bin) == 0) {
                /*
                 * Current rep has an empty bin, which is not OK.
                 * Mark this rep as dirty and break out of loop.
                 */
                clean = 0;
                break;
            }
            if(Tabulation_overflow(self->tab[rep])) {
                /*
                 * Current rep has invalid data because its Tabulation
                 * overflowed.
                 */
                clean = 0;
                break;
            }
        }
        /*
         * If current rep is clean, then increment reps; otherwise
         * free current rep, move in terminal rep, and reduce count of
         * good reps.
         */
        if(clean)
            ++rep;
        else {
            free(self->start[rep]);
            Tabulation_free(self->tab[rep]);
            Spectab_free(self->spectab[rep]);

            if(rep < nGoodReps - 1) {
                self->start[rep] = self->start[nGoodReps - 1];
                self->tab[rep] = self->tab[nGoodReps - 1];
                self->spectab[rep] = self->spectab[nGoodReps - 1];
            }
            --nGoodReps;
        }
    }
    self->nrep = nGoodReps;
#ifndef NDEBUG
    BootChr_sanityCheck(self, __FILE__, __LINE__);
#endif
    return self->nrep;
}

BootConf   *BootConf_new(BootChr * boot, double confidence) {
    BootConf   *bc = malloc(sizeof(BootConf));

    checkmem(bc, __FILE__, __LINE__);

    bc->confidence = confidence;
    bc->nrep = boot->nrep;
    bc->blockLength = boot->blockLength;
    bc->nBins = boot->nBins;

    bc->low = malloc(bc->nBins * sizeof(bc->low[0]));
    checkmem(bc->low, __FILE__, __LINE__);

    bc->high = malloc(bc->nBins * sizeof(bc->high[0]));
    checkmem(bc->high, __FILE__, __LINE__);

    bc->sep_cm = malloc(bc->nBins * sizeof(bc->sep_cm[0]));
    checkmem(bc->sep_cm, __FILE__, __LINE__);
    memset(bc->sep_cm, 0, bc->nBins * sizeof(bc->sep_cm[0]));

    assert(boot->spectab && boot->spectab[0]);

    // Dimension of site frequency spectrum
    bc->specDim = Spectab_dim(boot->spectab[0]);

    bc->loSpec = malloc(bc->specDim * sizeof(bc->loSpec[0]));
    checkmem(bc->loSpec, __FILE__, __LINE__);

    bc->hiSpec = malloc(bc->specDim * sizeof(bc->hiSpec[0]));
    checkmem(bc->hiSpec, __FILE__, __LINE__);

    int         i, rep, nvals;
    double      v[bc->nrep];

    // Confidence bounds on sigdsq
    for(i = 0; i < bc->nBins; ++i) {
        double      tmp1, nobs = 0, sigdsq;
        long unsigned tmp2;
        nvals = 0;

        for(rep = 0; rep < bc->nrep; ++rep) {
            sigdsq = Tabulation_sigdsq(boot->tab[rep], i, &tmp1, &tmp2);
            if(isfinite(sigdsq)) {
                assert(nvals < bc->nrep);
                v[nvals] = sigdsq;
                bc->sep_cm[i] += tmp1;
                nobs += tmp2;
                ++nvals;
            }
        }
        nobs /= nvals;
        bc->sep_cm[i] /= nvals;
        if(nvals < 10) {
            bc->low[i] = bc->high[i] = strtod("NAN", NULL);
        } else
            confidenceBounds(bc->low + i, bc->high + i, confidence, v, nvals);
    }

    // Confidence bounds on spectrum
    for(i = 0; i < bc->specDim; ++i) {
        long unsigned count;
        nvals = 0;

        for(rep = 0; rep < bc->nrep; ++rep) {
            count = Spectab_get(boot->spectab[rep], i);
            v[nvals++] = (double) count;
        }
        if(nvals < 10) {
            bc->low[i] = bc->high[i] = strtod("NAN", NULL);
        } else
            confidenceBounds(bc->loSpec + i, bc->hiSpec + i,
                             confidence, v, nvals);
    }

    return bc;
}

void BootConf_printHdr(const BootConf * self, FILE * ofp) {
    fprintf(ofp, "#%12s: %lg%% confidence bounds"
            " based on moving blocks bootstrap\n",
            "loLD, hiLD", 100.0 * self->confidence);
    fprintf(ofp, "#%12s: nrep=%ld blockLength=%ld nBins=%d\n",
            "Bootstrap parameters", self->nrep, self->blockLength, self->nBins);
}

double BootConf_lowBound(const BootConf * self, long bin) {
    return self->low[bin];
}

double BootConf_highBound(const BootConf * self, long bin) {
    return self->high[bin];
}

double BootConf_loSpecBound(const BootConf * self, long i) {
    return self->loSpec[i];
}

double BootConf_hiSpecBound(const BootConf * self, long i) {
    return self->hiSpec[i];
}

void BootConf_print(const BootConf * self, FILE * ofp) {
    int         i;

    BootConf_printHdr(self, ofp);
    fprintf(ofp, "%5s %10s %10s\n", "bin", "loLD", "hiLD");
    for(i = 0; i < self->nBins; ++i)
        fprintf(ofp, "%5d %10g %10g\n", i, self->low[i], self->high[i]);

    putc('\n', ofp);
    fprintf(ofp, "%5s %10s %10s\n", "i", "loSpec", "hiSpec");
    for(i = 0; i < self->specDim; ++i)
        fprintf(ofp, "%5d %10g %10g\n", i + 1, self->loSpec[i], self->hiSpec[i]);
}

void BootConf_free(BootConf * self) {
    free(self->low);
    free(self->high);
    free(self->sep_cm);
    free(self);
}

#endif
