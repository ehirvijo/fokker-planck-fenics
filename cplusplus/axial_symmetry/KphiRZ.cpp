
#include "KphiRZ.h"

#include <boost/math/special_functions/ellint_1.hpp>
using boost::math::ellint_1;

using namespace dolfin;

// A method to compute the value
void KphiRZ::eval(Array<double>& values, const Array<double>& x) const
{
  double mypi(3.14159265359);
  double k(4*r0*x[0]/((x[1]-z0)*(x[1]-z0)+(r0+x[0])*(r0+x[0])));
  
  if(fabs(k-1.)>1.e-12){
    values[0] = -ellint_1(sqrt(k))/(mypi*sqrt((x[1]-z0)*(x[1]-z0)+(r0+x[0])*(r0+x[0])));
  } else {
    values[0]=0.;
  }
  
  
}

