double        BranchTab_sum(const BranchTab *self);
int           BranchTab_normalize(BranchTab *self);
double        BranchTab_KLdiverg(const BranchTab *obs, const BranchTab *expt);

/// Return sum of values in BranchTab.
double BranchTab_sum(const BranchTab *self) {
    unsigned i;
    double s=0.0;

    for(i = 0; i < BT_DIM; ++i) {
        BTLink *el;
        for(el = self->tab[i]; el; el = el->next)
            s += el->value;
    }

    return s;
}

/// Divide all values by their sum. Return 0
/// on success, or 1 on failure.
int BranchTab_normalize(BranchTab *self) {
    unsigned i;
    double s = BranchTab_sum(self);

    if(s==0) 
        return 1;

    // divide by sum
    for(i = 0; i < BT_DIM; ++i) {
        BTLink *el;
        for(el = self->tab[i]; el; el = el->next)
            el->value /= s;
    }

    return 0;
}

/// Calculate KL divergence from two BranchTab objects, which
/// should be normalized before entering this function. Use
/// BranchTab_normalize to normalize. Function returns HUGE_VAL if there
/// are observed values without corresponding values in expt.
double BranchTab_KLdiverg(const BranchTab *obs, const BranchTab *expt) {
    assert(Dbl_near(1.0, BranchTab_sum(obs)));
	if(!Dbl_near(1.0, BranchTab_sum(expt))) {
		double s = BranchTab_sum(expt);
		double e = s-1.0;
		printf("%s:%s:%d: bad BranchTab sum=%lf err=%le\n",
			   __FILE__,__func__,__LINE__,
			   s, e);
		BranchTab_print(expt);
		fflush(stdout);
	}
    assert(Dbl_near(1.0, BranchTab_sum(expt)));

    int i;
    double kl=0.0;
    double p;  // observed frequency
    double q;  // frequency under model
    for(i=0; i < BT_DIM; ++i) {
        BTLink *o, *e;
        o = obs->tab[i];
        e = expt->tab[i];
        while(o && e) {
            if(o->key < e->key) {
                // e->value is q=0. This case blows up unless p==0.
                p = o->value;
                if(p != 0.0) {
                    kl = HUGE_VAL;
#if 0
                    fprintf(stderr,"%s:%s:%d: missing expt for ",
                            __FILE__,__func__,__LINE__);
                    printBits(sizeof(o->key), &o->key, stderr);
                    //exit(EXIT_FAILURE);
#endif
                    o = o->next;
                } else {
                    o = o->next;
                    continue;
                }
            }else if(o->key > e->key) {
                // o->value is p=0. For this case, note that
                // p*log(p/q) is the log of (p^p / q^p), which -> 1 as
                // p->0.  The log of 1 is 0, so 0 is the contribution
                // to kl.
                e = e->next;
                continue;
            }else {
                assert(o->key == e->key);
                p = o->value;
                q = e->value;
                e = e->next;
                o = o->next;
                kl += p*log(p/q);
            }
        }
        while(o) { // e->value is q=0: fail unless p=0
            p = o->value;
            if(p != 0.0) {
                kl = HUGE_VAL;
#if 0
                fprintf(stderr,"%s:%s:%d: missing expt for ",
                        __FILE__,__func__,__LINE__);
                printBits(sizeof(o->key), &o->key, stderr);
                exit(EXIT_FAILURE);
#endif
                o = o->next;
            } else {
                o = o->next;
                continue;
            }
        }
        // Any remaining cases have p=0, so contribution to kl is 0.
    }
    return kl;
}
