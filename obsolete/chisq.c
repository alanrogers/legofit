/// Bisect to find Chi-squared statistic that implies given value of
/// upper tail probability.
double getChiSqGoal(double df, double upTailProb) {
    double lo = 0.0, mid, hi = 3.0*df;
    double Qlo = gsl_sf_gamma_inc_Q(0.5*df, 0.5*lo);
    double Qhi;

    do{
        hi += 2.0*df;
        Qhi = gsl_sf_gamma_inc_Q(0.5*df, 0.5*hi);
    }while(Qhi > upTailProb);

    // Check input
    if( upTailProb > Qlo || upTailProb < Qhi ) {
        fprintf(stderr,"%s:%s:%d: initial interval doesn't enclose goal.\n",
                __FILE__,__func__,__LINE__);
        fprintf(stderr,"  lo=%lg Qlo=%lg hi=%lg Qhi=%lg goal=%lg\n",
                lo, Qlo, hi, Qhi, upTailProb);
        exit(EXIT_FAILURE);
    }

    // Bisect
    while(hi-lo > DBL_EPSILON*hi) {
        mid = lo + 0.5*(hi-lo);
        double Q = gsl_sf_gamma_inc_Q(0.5*df, 0.5*mid);
        if(Q < upTailProb)
            hi = mid;
        else
            lo = mid;
    }
    return hi;
}
