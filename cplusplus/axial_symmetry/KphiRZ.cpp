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

  // std::cout<<"r0 "<<r0<<" x0 "<<x[0]<<" x1 "<<x[1]<<" z0 "<<z0<<std::endl;
  // std::cout<<"k "<<k<<std::endl;
  
  if(fabs(k-1.)>1.e-12){
    values[0] = -ellint_1(sqrt(k))/(mypi*sqrt((x[1]-z0)*(x[1]-z0)+(r0+x[0])*(r0+x[0])));
    //    std::cout<<"elint "<<ellint_1(k)<<" value "<<values[0]<<std::endl;
  } else {
    values[0]=0.;
  }
  

}

