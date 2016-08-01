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
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

/// Contains the all data involved in a moving blocks bootstrap.
struct Boot {
    long        blockLength;    // number of SNPs per block
    long        nchr;           // number of chromosomes in data
    long        nrep;           // number of bootstrap replicates 
    int         npat;           // number of site patterns
    long       *nsnp;           // nsnp[i]: # snps on chr i
    long       *nblock;         // nblock[i]: blocks on chr i
    double     *patCount;       // patCount[i]: count of i'th site pattern
    long      **start;          // start[i][j] = start of j'th block in i'th rep
};

/** Contains the data for a bootstrap confidence interval. */
struct BootConf {
    long        nrep;           // repetitions
    long        blockLength;    // nucleotide positions per block
    double      confidence;     // size of confidence region
    double     *low, *high;     // confidence bounds
};

/// Constructor for class Boot.
Boot       *Boot_new(long nchr, long nsnp[nchr], long nrep, 
                     long blockLength, gsl_rng * rng) {
    long i, j;
    assert(blockLength > 0);
    if(nrep == 0)
        return NULL;

    for(i=0; i < nchr; ++i) {
        if(blockLength > nsnp[i]) {
            fprintf(stderr,
                    "ERR@%s:%d: in Boot_new, nsnp must be >blockLength.\n"
                    " Instead, nsnp[%ld]=%ld, blockLength=%ld.\n",
                    __FILE__, __LINE__, nsnp[i], blockLength);
            fprintf(stderr, " Use --blocksize argument"
                    " to reduce blockLength.\n");
            exit(1);
        }
    }

    // Block positions on chr i are uniform on [0,
    // nsnp[i]-blockLength+1).
    unsigned long endpos[nchr];
    for(i=0; i < nchr; ++i)
        endpos[i] = nsnp[i] - blockLength + 1;

    Boot       *boot = malloc(sizeof(Boot));
    CHECKMEM(boot);

    boot->nchr = nchr;
    boot->nsnp = memdup(nsnp, nchr*sizeof(boot->nsnp[i]));
    CHECKMEM(boot->nsnp);
    boot->nrep = nrep;
    boot->blockLength = blockLength;

    boot->nblock = LInt_div_round(nsnp, blockLength);
    boot->nBins = nBins;

    Boot_allocArrays(boot);

    for(i = 0; i < boot->nrep; ++i) {
        boot->tab[i] = Tabulation_new(windowcm, nBins);
        checkmem(boot->tab[i], __FILE__, __LINE__);

        boot->spectab[i] = Spectab_new(twoNsmp, folded);
        checkmem(boot->spectab[i], __FILE__, __LINE__);

        for(j = 0; j < boot->nblock; ++j)
            boot->start[i][j] = gsl_rng_uniform_int(rng, endpos);

        qsort(boot->start[i], (size_t) boot->nblock,
              sizeof(boot->start[0][0]), compareLongs);
    }

#ifndef NDEBUG
    Boot_sanityCheck(boot, __FILE__, __LINE__);
#endif
    return boot;
}

#if 0

#ifndef NDEBUG
void Boot_sanityCheck(const Boot * boot, const char *file, int line) {
    long        i, j;
    REQUIRE(boot != NULL, file, line);
    REQUIRE(boot->nchr > 0, file, line);
    REQUIRE(boot->nsnp != NULL, file, line);
    for(i=0; i < boot->nchr; ++i) {
        REQUIRE(boot->nsnp[i] > 0, file, line);
        REQUIRE(boot->nblock[i] > 0, file, line);
    }
    REQUIRE(boot->nrep > 0, file, line);
    REQUIRE(boot->blockLength > 0, file, line);
    REQUIRE(boot->blockLength < 100000, file, line);
    REQUIRE(boot->npat > 0, file, line);
    for(i=0; i < boot->npat; ++i) 
        REQUIRE(boot->patCount[i] >= 0.0, file, line);

    REQUIRE(boot->start != NULL, file, line);

    unsigned long endpos = boot->nsnp - boot->blockLength + 1;
    long        prev;

    for(i = 0; i < boot->nrep; ++i) {
        REQUIRE(boot->start[i] != NULL, file, line);
        for(j = 0; j < boot->nblock; ++j) {
            REQUIRE(boot->start[i][j] >= 0, file, line);
            REQUIRE(boot->start[i][j] < endpos, file, line);
            if(j > 0)
                REQUIRE(boot->start[i][j] >= prev, file, line);
            prev = boot->start[i][j];
        }
    }
}
#endif

/**
 * Allocate Boot's arrays. This code is used in several places, and I
 * have it here to ensure consistency.
 */
static void Boot_allocArrays(Boot * boot) {

    long        i;

    boot->start = calloc((unsigned long) boot->nrep, sizeof(boot->start[0]));
    checkmem(boot->start, __FILE__, __LINE__);

    boot->tab = calloc((unsigned long) boot->nrep, sizeof(boot->tab[0]));
    checkmem(boot->tab, __FILE__, __LINE__);

    boot->spectab = calloc((unsigned long) boot->nrep,
                           sizeof(boot->spectab[0]));
    checkmem(boot->spectab, __FILE__, __LINE__);

    for(i = 0; i < boot->nrep; ++i) {
        boot->start[i] = calloc((unsigned long) boot->nblock,
                                sizeof(boot->start[0][0]));
        checkmem(boot->start[i], __FILE__, __LINE__);
    }
}

/**
 * How many copies of snp are present in a given repetition (rep)?
 */
long Boot_multiplicity(const Boot * boot, long snp, long rep) {
    long        lndx, hndx, lowtarget;

    myassert(snp < boot->nsnp);

    /* lndx is index of first block containing snp */
    lowtarget = snp - boot->blockLength + 1;
    lndx = long_first_geq(lowtarget, boot->start[rep], boot->nblock);
    if(lndx == boot->nblock || boot->start[rep][lndx] > snp)
        return 0;

    myassert(snp >= boot->start[rep][lndx]);
    myassert(snp - boot->start[rep][lndx] < boot->blockLength);

    /* hndx is index of first block not containing snp */
    hndx = long_first_geq(snp + 1, boot->start[rep] + lndx,
                          boot->nblock - lndx);
    hndx += lndx;

    myassert(hndx == 0
             || boot->start[rep][hndx - 1] - snp < boot->blockLength);

    return hndx - lndx;
}

/*
 * Add one LD value to a Boot structure. On entry, boot points to the
 * Boot structure, Dsq and pqpq are the contributions to the numerator
 * and denominator of sigma_d^2, sep_cm is the separation in kilobases
 * between the two loci, and ndx1 and ndx2 are their positions within
 * the list of SNPs.
 */
void Boot_addLD(Boot * boot, double Dsq, double pqpq, double sep_cm,
                const SNP * snp1, const SNP * snp2) {
    for(register int rep = 0; rep < boot->nrep; ++rep) {
        register unsigned wgt = SNP_multiplicity(snp1, rep)
            * SNP_multiplicity(snp2, rep);

        if(wgt == 0)
            continue;
        Tabulation_record(boot->tab[rep], Dsq, pqpq, sep_cm, wgt);
    }
}

/*
 * Allocate and initialize a duplicate of a Boot structure. Deallocate
 * with Boot_free.
 */
Boot       *Boot_dup(const Boot * old) {
    myassert(old != NULL);
    myassert(old->nblock > 0);
    Boot       *new = memdup(old, sizeof(*old));

    myassert(new->nblock > 0);

    new->start = malloc(new->nrep * sizeof(new->start[0]));
    checkmem(new->start, __FILE__, __LINE__);

    new->tab = malloc(new->nrep * sizeof(new->tab[0]));
    checkmem(new->tab, __FILE__, __LINE__);

    new->spectab = malloc(new->nrep * sizeof(new->spectab[0]));
    checkmem(new->spectab, __FILE__, __LINE__);

    myassert(sizeof(new->start[0][0]) > 0);

    for(int i = 0; i < new->nrep; ++i) {
        new->tab[i] = Tabulation_dup(old->tab[i]);
        new->spectab[i] = Spectab_dup(old->spectab[i]);
        new->start[i] = memdup(old->start[i],
                               new->nblock * sizeof(new->start[0][0]));
    }

#ifndef NDEBUG
    Boot_sanityCheck(new, __FILE__, __LINE__);
#endif
    return new;
}

/* Return number of bootstrap repetitions */
long Boot_nrep(const Boot * boot) {
    myassert(boot);

    return boot->nrep;
}

/* Return number of bins */
int Boot_nBins(const Boot * boot) {
    myassert(boot);

    return boot->nBins;
}

/** Return number of blocks */
long Boot_nblock(const Boot * boot) {
    myassert(boot);

    return boot->nblock;
}

/** Return number of SNPs */
long Boot_nsnp(const Boot * boot) {
    myassert(boot);

    return boot->nsnp;
}

/* deallocate bootstrap */
void Boot_free(Boot * boot) {
    long        i;

#ifndef NDEBUG
    Boot_sanityCheck(boot, __FILE__, __LINE__);
#endif

    for(i = 0; i < boot->nrep; ++i) {
        free(boot->start[i]);
        Tabulation_free(boot->tab[i]);
        Spectab_free(boot->spectab[i]);
        boot->start[i] = NULL;
        boot->tab[i] = NULL;
    }
    free(boot->start);
    free(boot->tab);
    free(boot->spectab);
    boot->start = NULL;
    boot->tab = NULL;
    free(boot);
}

/**
 * Interpolate in order to approximate the value v[p*(len-1)].
 * Return NaN if len==0.
 */
double interpolate(double p, double *v, long len) {
    if(len == 0)
        return strtod("NAN", 0);
    long        i, j;
    double      w;
    double      goal = p * (len - 1);

    i = floor(goal);
    j = ceil(goal);

    myassert(i >= 0);
    myassert(j < len);

    if(i == j)                  /* no interpolation needed */
        return v[i];
    w = goal - i;
    return (1.0 - w) * v[i] + w * v[j];
}

/* Return 1 if the two Boot structs are idential; zero otherwise */
int Boot_equals(const Boot * x, const Boot * y) {
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
 * Add the contents of Boot structure y to Boot structure x.
 * On return x will summarize all the LD values that were originally
 * contained in the two original Boot structures.
 */
void Boot_plus_equals(Boot * x, const Boot * y) {
    int         rep;

#ifndef NDEBUG
    Boot_sanityCheck(x, __FILE__, __LINE__);
    Boot_sanityCheck(y, __FILE__, __LINE__);
#endif
    if(x == NULL)
        die("Boot_plus_equals: destination is NULL", __FILE__, __LINE__);
    if(y == NULL)
        return;

    if(x->nsnp != y->nsnp
       || x->nrep != y->nrep
       || x->blockLength != y->blockLength
       || x->nblock != y->nblock || x->nBins != y->nBins)
        die("Boot_plus_equals: unconformable arguments", __FILE__, __LINE__);

    for(rep = 0; rep < x->nrep; ++rep) {
        long        block;

        for(block = 0; block < x->nblock; ++block) {
            if(x->start[rep][block] != y->start[rep][block])
                die("Boot_plus_equals: incompatible block starts",
                    __FILE__, __LINE__);
        }
        Tabulation_plus_equals(x->tab[rep], y->tab[rep]);
        Spectab_plus_equals(x->spectab[rep], y->spectab[rep]);
    }
#ifndef NDEBUG
    Boot_sanityCheck(x, __FILE__, __LINE__);
#endif
    return;
}

/*
 * Write a text file representing a Boot structure. To restore the
 * Boot structure from this file, use Boot_restore.
 */
void Boot_dump(const Boot * boot, FILE * ofp) {
    long        rep, block;

    if(fprintf(ofp, " %ld %ld %ld %ld %d\n",
               boot->nsnp, boot->nrep, boot->blockLength,
               boot->nblock, boot->nBins) < 0)
        eprintf("fprintf", __FILE__, __LINE__);

    for(rep = 0; rep < boot->nrep; ++rep) {
        for(block = 0; block < boot->nblock; ++block) {
            if(fprintf(ofp, " %ld", boot->start[rep][block]) < 0)
                eprintf("fprintf", __FILE__, __LINE__);
        }
        putc('\n', ofp);
        Tabulation_dump(boot->tab[rep], ofp);
        Spectab_dump(boot->spectab[rep], ofp);
    }
}

/*
 * Create a Boot structure by reading a file produced by
 * Boot_dump.
 */
Boot       *Boot_restore(FILE * ifp) {

    long        nsnp, nrep, blockLength, nblock, rep, block;
    int         nBins, rval;
    Boot       *boot = malloc(sizeof(Boot));

    checkmem(boot, __FILE__, __LINE__);

    rval = fscanf(ifp, "%ld %ld %ld %ld %d",
                  &nsnp, &nrep, &blockLength, &nblock, &nBins);
    if(rval != 5)
        eprintf("ERR@%s:%d: fscanf returned %d instead of 5\n",
                __FILE__, __LINE__, rval);

    boot->nsnp = nsnp;
    boot->nrep = nrep;
    boot->blockLength = blockLength;
    boot->nblock = (long) round(nsnp / ((double) blockLength));
    boot->nBins = nBins;

    if(nblock != boot->nblock) {
        fprintf(stderr, "from file:"
                " nsnp=%ld nrep=%ld blockLength=%ld nblock=%ld nBins=%d\n",
                nsnp, nrep, blockLength, nblock, nBins);
        fprintf(stderr, "in boot  :"
                " nsnp=%ld nrep=%ld blockLength=%ld nblock=%ld nBins=%d\n",
                boot->nsnp, boot->nrep, boot->blockLength,
                boot->nblock, boot->nBins);
        eprintf("ERR@%s:%d: nblock=%ld but boot->nblock=%ld\n",
                __FILE__, __LINE__, nblock, boot->nblock);
    }

    Boot_allocArrays(boot);

    for(rep = 0; rep < boot->nrep; ++rep) {
        for(block = 0; block < boot->nblock; ++block) {
            rval = fscanf(ifp, " %ld", boot->start[rep] + block);
            if(rval != 1)
                eprintf("ERR@%s:%d: fscanf returned %d instead of 1\n",
                        __FILE__, __LINE__, rval);
        }
        boot->tab[rep] = Tabulation_restore(ifp);
        boot->spectab[rep] = Spectab_restore(ifp);
        if(!Tabulation_isfinite(boot->tab[rep])) {
            fprintf(stderr,
                    "Warning@%s:%d: bootstrap replicate %ld is non-finite\n",
                    __FILE__, __LINE__, rep);
        }
    }
#ifndef NDEBUG
    Boot_sanityCheck(boot, __FILE__, __LINE__);
#endif
    return boot;
}

/**
 * Fill arrays sigdsq, cm, and nobs with values for bootstrap
 * repetition "rep". If nobs==NULL, nothing is stored there.
 */
void Boot_get_rep(Boot * boot, DblArray * sigdsq, DblArray * rsq,
                  DblArray * cm, ULIntArray * nobs,
                  ULIntArray * spectrum, int rep) {
    myassert(boot);
    myassert(sigdsq);
    myassert(cm);
    myassert(rep < boot->nrep);
    myassert(rep >= 0);
    Tabulation_report(boot->tab[rep], cm, nobs, sigdsq, rsq);
    Spectab_report(boot->spectab[rep], spectrum);
}

/*
 * Return the raw counts for a single rep and bin. Counts for
 * numerator, denominator, and sep_cm are returned in the
 * corresponding arguments. The function itself returns nobs.
 */
long unsigned Boot_rawCounts(const Boot * boot, int rep, int bin,
                             double *numerator, double *denominator,
                             double *sumRsq, double *sep_cm) {

    assert(rep < boot->nrep);
    assert(bin < boot->nBins);

    return Tabulation_rawCounts(boot->tab[rep], bin, numerator,
                                denominator, sumRsq, sep_cm);
}

/**
 * Remove replicates in which some bins have zero observations.
 * These would generate NaN values in the calculation of sigdsq,
 * which are not handled by the minimizer. Return revised number of
 * bootstrap replicates.
 */
long Boot_purge(Boot * boot) {
    myassert(boot);

    long        rep, nGoodReps = boot->nrep;
    int         bin;

    rep = 0;
    while(rep < nGoodReps) {
        int         clean = 1;

        for(bin = 0; bin < boot->nBins; ++bin) {
            if(Tabulation_nObs(boot->tab[rep], bin) == 0) {
                /*
                 * Current rep has an empty bin, which is not OK.
                 * Mark this rep as dirty and break out of loop.
                 */
                clean = 0;
                break;
            }
            if(Tabulation_overflow(boot->tab[rep])) {
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
            free(boot->start[rep]);
            Tabulation_free(boot->tab[rep]);
            Spectab_free(boot->spectab[rep]);

            if(rep < nGoodReps - 1) {
                boot->start[rep] = boot->start[nGoodReps - 1];
                boot->tab[rep] = boot->tab[nGoodReps - 1];
                boot->spectab[rep] = boot->spectab[nGoodReps - 1];
            }
            --nGoodReps;
        }
    }
    boot->nrep = nGoodReps;
#ifndef NDEBUG
    Boot_sanityCheck(boot, __FILE__, __LINE__);
#endif
    return boot->nrep;
}

/** Print a Boot object */
void Boot_print(const Boot * boot, FILE * ofp) {
    long        rep, block;

    fprintf(ofp,
            "Boot_print: nsnp=%ld nrep=%ld blockLength=%ld nblock=%ld\n",
            boot->nsnp, boot->nrep, boot->blockLength, boot->nblock);

    fprintf(ofp, "Block starts:\n");
    for(rep = 0; rep < boot->nrep; ++rep) {
        fprintf(ofp, "  rep %ld:", rep);
        for(block = 0; block < boot->nblock; ++block)
            fprintf(ofp, " %ld", boot->start[rep][block]);
        putc('\n', ofp);
        // Tabulation_print(boot->tab[rep], ofp);
        // Spectab_print(boot->spectab[rep], ofp);
    }
}

#ifndef NDEBUG
/* For debugging Boot_multiplicity */
unsigned Boot_multiplicity_slow(Boot * boot, long snp, long rep) {
    unsigned    i, n = 0;

    for(i = 0; i < boot->nblock; ++i) {
        long        distance = snp - boot->start[rep][i];

        if(distance < 0)
            break;
        if(distance < boot->blockLength)
            ++n;
    }
    return n;
}
#endif

BootConf   *BootConf_new(Boot * boot, double confidence) {
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
                myassert(nvals < bc->nrep);
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

void BootConf_printHdr(const BootConf * bc, FILE * ofp) {
    fprintf(ofp, "#%12s: %lg%% confidence bounds"
            " based on moving blocks bootstrap\n",
            "loLD, hiLD", 100.0 * bc->confidence);
    fprintf(ofp, "#%12s: nrep=%ld blockLength=%ld nBins=%d\n",
            "Bootstrap parameters", bc->nrep, bc->blockLength, bc->nBins);
}

double BootConf_lowBound(const BootConf * bc, long bin) {
    return bc->low[bin];
}

double BootConf_highBound(const BootConf * bc, long bin) {
    return bc->high[bin];
}

double BootConf_loSpecBound(const BootConf * bc, long i) {
    return bc->loSpec[i];
}

double BootConf_hiSpecBound(const BootConf * bc, long i) {
    return bc->hiSpec[i];
}

void BootConf_print(const BootConf * bc, FILE * ofp) {
    int         i;

    BootConf_printHdr(bc, ofp);
    fprintf(ofp, "%5s %10s %10s\n", "bin", "loLD", "hiLD");
    for(i = 0; i < bc->nBins; ++i)
        fprintf(ofp, "%5d %10g %10g\n", i, bc->low[i], bc->high[i]);

    putc('\n', ofp);
    fprintf(ofp, "%5s %10s %10s\n", "i", "loSpec", "hiSpec");
    for(i = 0; i < bc->specDim; ++i)
        fprintf(ofp, "%5d %10g %10g\n", i + 1, bc->loSpec[i], bc->hiSpec[i]);
}

void BootConf_free(BootConf * bc) {
    free(bc->low);
    free(bc->high);
    free(bc->sep_cm);
    free(bc);
}

#endif
