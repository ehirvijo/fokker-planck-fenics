// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#include "GSource.h"

using namespace dolfin;

// Evaluate the value of the source at requested point
// Remember that the coordinates are x[0] = r, x[1] = z
void GSource::eval(Array<double>& values, const Array<double>& x) const
{

  values[0] = -g0*fabs(x[1]);  // note: symmetric around v||=0
}

// A method to compute the coefficents at this time
void GSource::compute_coeffs()
{
  g0 = gmax*(1.0 - exp(-gamma_g*t));
}

