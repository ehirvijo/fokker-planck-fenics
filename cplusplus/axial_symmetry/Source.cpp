// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#include "Source.h"

using namespace dolfin;

// Evaluate the value of the source at requested point
// Remember that the coordinates are x[0] = r, x[1] = z
void Source::eval(Array<double>& values, const Array<double>& x) const
{
  values[0] = 0.0*exp(-(x[0]*x[0]+x[1]*x[1])/0.01);
}

