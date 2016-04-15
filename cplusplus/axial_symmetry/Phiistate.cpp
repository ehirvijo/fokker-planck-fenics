
#include "Phiistate.h"

using namespace dolfin;

// A method to compute the value
void Phiistate::eval(Array<double>& values, const Array<double>& x) const
{
  double mypi(3.14159265359);
  double rx2 = x[0]*x[0]+x[1]*x[1];

  // Limit[Erf[x]/x,x->0]=2/sqrt(mypi)
  if (sqrt(rx2)<1.e-6) { 
    values[0] = -phii0/2.0;
  } else {
    values[0] = -phii0*sqrt(mypi)*erf(sqrt(rx2))/(4.0*sqrt(rx2));
  }
  
}

// A method to compute the coefficents at this time
void Phiistate::compute_coeffs() 
{
  phii0 = phimax*((1.0 + exp(-gamma_phi*t))/2.0);
}
