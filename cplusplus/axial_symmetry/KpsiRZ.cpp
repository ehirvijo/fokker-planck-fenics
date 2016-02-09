// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#include "KpsiRZ.h"

#include <boost/math/special_functions/ellint_2.hpp>
using boost::math::ellint_2;

using namespace dolfin;

// A method for setting the parameters
void KphiRZ::set(const Array<double>& x)
{
  r0=x[0];
  z0=x[1];
}
 
// A method to compute the value
void KphiRZ::eval(Array<double>& values, const Array<double>& x) const
{
  double mypi(3.14159265359);
  double k(4*r0*x[0]/((x[1]-z0)*(x[1]-z0)+(r0+x[0])*(r0+x[0])));

  values[0] = -ellint_2(sqrt(k))*sqrt((x[1]-z0)*(x[1]-z0)+(r0+x[0])*(r0+x[0]))/(2*mypi);

}

