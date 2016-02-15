// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#include "KpsiRZ.h"

#include <boost/math/special_functions/ellint_2.hpp>
using boost::math::ellint_2;

using namespace dolfin;
 
// A method to compute the value
void KpsiRZ::eval(Array<double>& values, const Array<double>& x) const
{
  double mypi(3.14159265359);
  double k(4*r0*x[0]/((x[1]-z0)*(x[1]-z0)+(r0+x[0])*(r0+x[0])));

  // construct a check for the problematic case when 
  if (fabs((r0+x[0])*(r0+x[0])+(x[1]-z0)*(x[1]-z0))<=1.e-15 ) {
    values[0] = 0.;
  } else {
    //std::cout<<"r0="<<r0<<" z0="<<z0<<" r'="<<x[0]<<" z'="<<x[1]<<" k: "<<k<<std::endl;
    values[0] = -ellint_2(sqrt(k))*sqrt((x[1]-z0)*(x[1]-z0)+(r0+x[0])*(r0+x[0]))/(2*mypi);
  }
}

