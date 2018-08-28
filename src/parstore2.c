struct Param {
    char *name;
    double *value;
    double min, max;
};


struct ParStore2 {
    int         nFixed, nFree, nGaussian, nConstrained; // num pars
    double      loFree[MAXPAR]; // lower bounds
    double      hiFree[MAXPAR]; // upper bounds
    char       *nameFixed[MAXPAR];  // Parameter names
    char       *nameFree[MAXPAR];   // Parameter names
    char       *nameGaussian[MAXPAR];    // Parameter names
    char       *nameConstrained[MAXPAR]; // Parameter names
    ParKeyVal  *pkv;           // linked list of name/ptr pairs
    double      fixedVal[MAXPAR];   // parameter values
    double      freeVal[MAXPAR];    // parameter values
    double      gaussianVal[MAXPAR]; // parameter values
    double      constrainedVal[MAXPAR]; // parameter values
    double      mean[MAXPAR];        // Gaussian means
    double      sd[MAXPAR];          // Gaussian standard deviations
    te_expr    *constr[MAXPAR];      // controls constrainedVal entries
    te_variable te_pars[MAXPAR];
    char       *formulas[MAXPAR];    // formulas of constrained vars
};
