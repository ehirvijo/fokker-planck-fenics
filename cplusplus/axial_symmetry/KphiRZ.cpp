// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#include "KphiRZ.h"

#include <boost/math/special_functions/ellint_1.hpp>
using boost::math::ellint_1;

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
  values[0] = -ellint_1(k)/(mypi*sqrt((x[1]-z0)*(x[1]-z0)+(r0+x[0])*(r0+x[0])));
}


// default constructor for the type KphiRZ
KphiRZ::KphiRZ(): Expression(), r0(0.0),z0(0.0) {}

