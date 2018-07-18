/**
@file bma.c
@page bma
@author Alan R. Rogers
@brief Bootstrap model averaging

# `bma`: bootstrap model averaging

Bootstrap model averaging was proposed by Buckland et al (Biometrics,
53(2):603-618). It can be used with weights provided by any method of
model selection. Model selection is applied to the real data and also
to a set of bootstrap replicates. The weight, \f$w_i\f$ of the i'th
model is the fraction of these data sets for which the i'th model
wins.

The model-averaged estimator of a parameter, \f$\theta\f$, is the
average across models, weighted by \f$w_i\f$, of the model-specific
estimates of \f$\theta\f$. Models in which \f$\theta\f$ does not
appear are omitted from this weighted average.


@copyright Copyright (c) 2018, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/
