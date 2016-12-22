//  Pseudorandom numbers from a truncated Gaussian distribution.
//
// This code was translated by Alan R. Rogers into C, based on the C++
// original written by G. Dolle and V Mazet, which is available at
// http://miv.u-strasbg.fr/mazet/rtnorm/rtnormCpp.zip
//
//  This implements an extension of Chopin's algorithm detailed in
//  N. Chopin, "Fast simulation of truncated Gaussian distributions",
//  Stat Comput (2011) 21:275-288
//
//  Copyright (C) 2012 Guillaume Dollé, Vincent Mazet
//  (LSIIT, CNRS/Université de Strasbourg)
//  Version 2012-07-04, Contact: vincent.mazet@unistra.fr
//
//  06/07/2012:
//  - first launch of rtnorm.cpp
//
//  Licence: GNU General Public License Version 2
//  This program is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2 of the License, or (at your
//  option) any later version. This program is distributed in the hope that
//  it will be useful, but WITHOUT ANY WARRANTY; without even the implied
//  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details. You should have received a
//  copy of the GNU General Public License along with this program; if not,
//  see http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
//
//  Depends: LibGSL
//  OS: Unix based system


#ifndef DTNORM_H
#define DTNORM_H

#include <gsl/gsl_rng.h>

// Pseudorandom numbers from a truncated Gaussian distribution The
// Gaussian has parameters mu (default 0) and sigma (default 1) and is
// truncated on the interval [a,b].  Returns the random variate x.
double dtnorm(const double mu, const double sigma, double a, double b,
              gsl_rng *rng);

#endif
